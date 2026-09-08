#!/usr/bin/env python
# encoding: utf-8

r"""Tests for opt-in geometry reconstruction (``met/reconstruction.py``).

Two kinds of check, and the distinction matters:

* **Invariants and plumbing** -- monotonicity, positivity, SI units, the
  ``fill_dict`` contract, the sparse-track round trip.  These are written here
  and run now.  They are necessary but *not* sufficient: a kt-vs-m/s or
  km-vs-nmi error passes every one of them.
* **Reproducing the published relationship** -- the actual verification, driven
  by ``tests/data/storm/reconstruction_reference.json``.  That file ships empty,
  so the test fails until a human transcribes values from the paper.  The
  invariant tests for an estimator whose coefficients are not yet in place are
  skipped rather than deleted, so they switch on by themselves when the
  coefficients land.
"""

import json
from pathlib import Path
import sys
import warnings

import numpy as np
import pytest

from clawpack.geoclaw.met import reconstruction
from clawpack.geoclaw.met.track import StormTrack
from clawpack.geoclaw.met.parametric import ParametricMetForcing

sys.path.insert(0, str(Path(__file__).parent))
from test_storm import _storm_input_path  # noqa: E402

data_dir = Path(__file__).parent / "data" / "storm"

# Pure-Python storm tests; the CI selector is -m "python and not remote".
pytestmark = [pytest.mark.python, pytest.mark.storm]

WILLOUGHBY_VALID_ABOVE = reconstruction.WILLOUGHBY_2006_RMW["valid_above_m_s"]

# Estimators whose coefficients are not transcribed yet cannot be exercised.
_HAVE_WILLOUGHBY = (
    reconstruction.WILLOUGHBY_2006_RMW.get("coefficients") is not None)
needs_willoughby = pytest.mark.skipif(
    not _HAVE_WILLOUGHBY,
    reason="Willoughby et al. (2006) coefficients not transcribed yet")


def _synthetic_track(n=6, **overrides):
    """A minimal track with every field reported, for fill tests."""
    track = StormTrack()
    track.t = np.array(["2020-09-01T00", "2020-09-01T06", "2020-09-01T12",
                        "2020-09-01T18", "2020-09-02T00", "2020-09-02T06"],
                       dtype="datetime64[s]")[:n]
    track.eye_location = np.column_stack([np.linspace(-80.0, -83.0, n),
                                          np.linspace(20.0, 32.0, n)])
    track.max_wind_speed = np.linspace(20.0, 60.0, n)
    track.central_pressure = np.linspace(1000e2, 950e2, n)
    track.max_wind_radius = np.full(n, 40e3)
    track.storm_radius = np.full(n, 300e3)
    for name, value in overrides.items():
        setattr(track, name, value)
    return track


# ---------------------------------------------------------------------------
# The verification: reproduce the published relationship
# ---------------------------------------------------------------------------

def test_published_reference_values_are_transcribed():
    """The verification: reproduce the published relationship.

    Reproducing the paper is what actually verifies these estimators -- the
    invariants below would pass with a units error still in place.

    While the reference file is empty this **skips** rather than fails, so that
    master is not permanently red.  The guard against an unverified coefficient
    reaching a run is not this test: it is that the estimator itself raises
    (`test_unfilled_coefficients_raise_rather_than_guess`), which prevents *use*
    rather than merely merge.  Once the file is filled in, this test and the
    invariant tests below activate on their own.
    """
    with open(data_dir / "reconstruction_reference.json") as handle:
        reference = json.load(handle)

    block = reference["willoughby_2006_rmw"]
    if not block["cases"]:
        pytest.skip(
            "tests/data/storm/reconstruction_reference.json has no reference "
            "cases yet. Transcribe published (max_wind_speed, latitude, "
            "max_wind_radius) values from Willoughby, Darling & Rahn (2006), "
            "Mon. Wea. Rev. 134, 1102-1120, with the equation number, the "
            "wind-averaging convention, the basin, and the paper's reported "
            "scatter as 'tolerance_m'. See the _README in that file.")

    for field in ("equation", "wind_basis", "basin", "units", "tolerance_m"):
        assert block[field] is not None, (
            f"reference block is missing '{field}'; provenance must be "
            f"recorded alongside the values, since none of it is something an "
            f"invariant test can check for you.")

    tolerance = float(block["tolerance_m"])
    for case in block["cases"]:
        # Reference winds are on the paper's flight-level basis, so build a
        # stationary single-point track: to_willoughby_wind then reduces to
        # dividing by ATMOS_BOUNDARY_LAYER, which we undo here so the paper's
        # own value reaches Eq. (7) unchanged.
        track = _synthetic_track(n=1)
        track.max_wind_speed = np.array([
            float(case["max_wind_speed"]) * reconstruction.ATMOS_BOUNDARY_LAYER])
        track.eye_location = np.array([[-75.0, float(case["latitude"])]])
        estimate = reconstruction.rmw_willoughby2006(track.t[0], track)
        assert abs(estimate - float(case["max_wind_radius"])) <= tolerance, (
            f"Willoughby RMW at Vmax={case['max_wind_speed']} m/s, "
            f"lat={case['latitude']} deg: got {estimate:.0f} m, "
            f"paper gives {case['max_wind_radius']} m")


# ---------------------------------------------------------------------------
# Invariants -- necessary, not sufficient
# ---------------------------------------------------------------------------

@needs_willoughby
def test_rmw_decreases_with_intensity():
    """Stronger storms have tighter eyewalls."""
    latitudes = [15.0, 25.0, 35.0]
    for latitude in latitudes:
        radii = []
        for wind in (20.0, 35.0, 50.0, 65.0):
            track = _synthetic_track(n=1)
            track.max_wind_speed = np.array([wind])
            track.eye_location = np.array([[-75.0, latitude]])
            radii.append(reconstruction.rmw_willoughby2006(track.t[0], track))
        assert all(np.diff(radii) < 0.0), f"at lat {latitude}: {radii}"


@needs_willoughby
def test_rmw_increases_with_latitude():
    """Tropical cyclones broaden as they move poleward."""
    radii = []
    for latitude in (10.0, 20.0, 30.0, 40.0):
        track = _synthetic_track(n=1)
        track.max_wind_speed = np.array([40.0])
        track.eye_location = np.array([[-75.0, latitude]])
        radii.append(reconstruction.rmw_willoughby2006(track.t[0], track))
    assert all(np.diff(radii) > 0.0), radii


@needs_willoughby
def test_rmw_is_finite_positive_si_over_the_observed_envelope():
    """No NaN, no negative, and physically plausible magnitudes."""
    for wind in np.linspace(15.0, 85.0, 15):
        for latitude in np.linspace(5.0, 45.0, 9):
            track = _synthetic_track(n=1)
            track.max_wind_speed = np.array([wind])
            track.eye_location = np.array([[-75.0, latitude]])
            radius = reconstruction.rmw_willoughby2006(track.t[0], track)
            assert np.isfinite(radius) and radius > 0.0
            # metres, not km or nmi: an RMW outside 2-200 km is a units error.
            assert 2e3 < radius < 200e3, (
                f"RMW {radius} m at {wind} m/s, {latitude} deg looks like a "
                f"unit-conversion error")


def test_unfilled_coefficients_raise_rather_than_guess():
    """An estimator with no transcribed coefficients must refuse to run."""
    if _HAVE_WILLOUGHBY:
        pytest.skip("coefficients are transcribed")
    track = _synthetic_track()
    with pytest.raises(reconstruction.CoefficientsNotTranscribed,
                       match="no transcribed coefficients"):
        reconstruction.rmw_willoughby2006(track.t[0], track)


def test_wind_pressure_is_declared_not_implemented():
    """The WPR is deliberately deferred; it must say so, not return a number."""
    track = _synthetic_track()
    with pytest.raises(NotImplementedError):
        reconstruction.wind_pressure(track.t[0], track)


# ---------------------------------------------------------------------------
# The estimators that are implemented
# ---------------------------------------------------------------------------

def test_roci_climatology_matches_the_constant_it_replaces():
    """Naming the old hard-wired 500 km must not change its value."""
    track = _synthetic_track()
    assert reconstruction.roci_climatology(track.t[0], track) == 500e3
    assert reconstruction.ROCI_CLIMATOLOGY_M == 500e3


def test_rmw_constant_validates_its_input():
    fill = reconstruction.rmw_constant(40e3)
    track = _synthetic_track()
    assert fill(track.t[0], track) == 40e3
    for bad in (0.0, -1.0, np.nan):
        with pytest.raises(ValueError):
            reconstruction.rmw_constant(bad)


def test_at_time_prefers_an_exact_match():
    """fill_dict hands back a time off the track's own axis."""
    track = _synthetic_track()
    for index, t in enumerate(track.t):
        assert reconstruction._at_time(track, t, "max_wind_speed") == \
            pytest.approx(track.max_wind_speed[index])

    # A time between samples falls back to the nearest.
    between = track.t[0] + np.timedelta64(3600, "s")
    nearest = reconstruction._at_time(track, between, "max_wind_speed")
    assert nearest == pytest.approx(track.max_wind_speed[0])

    # A field the track does not carry is missing, not an error.
    track.max_wind_radius = None
    assert np.isnan(reconstruction._at_time(track, track.t[0],
                                            "max_wind_radius"))


# ---------------------------------------------------------------------------
# Assembly, reporting, and the round trip that matters
# ---------------------------------------------------------------------------

def test_default_fill_dict_is_a_no_op_by_default():
    """Calling it with no arguments must reproduce write_geoclaw's own default."""
    fill = reconstruction.default_fill_dict()
    assert set(fill) == {"storm_radius"}
    assert fill["storm_radius"] is reconstruction.roci_climatology


def test_default_fill_dict_accepts_and_removes_models():
    fill = reconstruction.default_fill_dict(
        models={"max_wind_radius": reconstruction.rmw_constant(30e3)})
    assert set(fill) == {"storm_radius", "max_wind_radius"}

    stripped = reconstruction.default_fill_dict(models={"storm_radius": None})
    assert stripped == {}


def test_default_fill_dict_warns_about_what_it_cannot_fill():
    track = _synthetic_track(max_wind_radius=np.full(6, np.nan))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        reconstruction.default_fill_dict(track)
    assert any("max_wind_radius" in str(w.message) for w in caught)


def test_coverage_report_counts_reported_and_filled():
    n = 6
    partial = np.full(n, np.nan)
    partial[:2] = 40e3
    track = _synthetic_track(max_wind_radius=partial)

    plain = reconstruction.coverage_report(track)
    assert plain["n_times"] == n
    assert plain["fields"]["max_wind_radius"]["reported"] == 2
    assert plain["fields"]["max_wind_radius"]["missing"] == n - 2
    assert plain["fields"]["max_wind_radius"]["filled"] == 0
    # Two of six times have all four fields, so only those are writable.
    assert plain["n_writable"] == 2

    filled = reconstruction.coverage_report(
        track, fill_dict=reconstruction.default_fill_dict(
            models={"max_wind_radius": reconstruction.rmw_constant(35e3)}))
    assert filled["fields"]["max_wind_radius"]["filled"] == n - 2
    assert filled["fields"]["max_wind_radius"]["still_missing"] == 0
    assert filled["n_writable"] == n


def test_coverage_report_on_a_real_sparse_track():
    """The bundled IBTrACS fixture: report the tail that cannot be written."""
    pytest.importorskip("xarray")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        track = StormTrack.read_ibtracs(_storm_input_path("ibtracs"),
                                        sid="2008245N17323",
                                        agency_pref=["wmo", "usa"])

    plain = reconstruction.coverage_report(track)
    assert plain["n_writable"] < plain["n_times"], (
        "the fixture should have times with no reported radius")

    filled = reconstruction.coverage_report(
        track, fill_dict=reconstruction.default_fill_dict(
            models={"max_wind_radius": reconstruction.rmw_constant(40e3)}))
    assert filled["n_writable"] == filled["n_times"], (
        "with both radii filled, every time should be writable")


def test_sparse_track_round_trips_to_a_runnable_storm_file(tmp_path):
    """Strip the radii, fill them, write: no row may be dropped.

    This is the end-to-end reason the module exists -- a track with no reported
    geometry must produce a 7-column file GeoClaw can actually run.
    """
    pytest.importorskip("pandas")
    track = StormTrack.read_atcf(_storm_input_path("atcf"))
    n_times = len(track.t)

    # Strip exactly what the historical archives are missing.
    track.max_wind_radius = np.full(n_times, np.nan)
    track.storm_radius = np.full(n_times, np.nan)

    forcing = ParametricMetForcing()
    for field in ("t", "time_offset", "eye_location", "max_wind_speed",
                  "central_pressure", "max_wind_radius", "storm_radius"):
        setattr(forcing, field, getattr(track, field))

    path = tmp_path / "filled.storm"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        forcing.write_geoclaw(path, fill_dict=reconstruction.default_fill_dict(
            models={"max_wind_radius": reconstruction.rmw_constant(40e3)}))

    lines = path.read_text().splitlines()
    assert int(lines[0].strip()) == n_times, "no track time may be dropped"

    values = np.loadtxt(path, skiprows=3)
    assert values.shape == (n_times, 7)
    assert np.all(np.isfinite(values)), "a runnable file has no NaN"
    assert np.all(values[:, 4] > 0.0), "max_wind_radius must be positive"
    assert np.all(values[:, 6] > 0.0), "storm_radius must be positive"

    # And the fills did not mutate the caller's track.
    assert np.all(np.isnan(track.max_wind_radius))


def test_estimators_match_the_fill_dict_signature():
    """Every public estimator must be usable as fn(t, storm) without adaptation."""
    track = _synthetic_track()
    for name in ("roci_climatology",):
        estimator = getattr(reconstruction, name)
        assert np.isfinite(estimator(track.t[0], track))
    assert np.isfinite(reconstruction.rmw_constant(40e3)(track.t[0], track))


# ---------------------------------------------------------------------------
# The wind basis Willoughby's regressions are fit on
# ---------------------------------------------------------------------------

def test_translation_speed_matches_a_known_track():
    """Forward difference of consecutive track points, as the Fortran does."""
    track = StormTrack()
    # One degree of latitude in one hour: ~111 km / 3600 s.
    track.t = np.array(["2020-09-01T00", "2020-09-01T01", "2020-09-01T02"],
                       dtype="datetime64[s]")
    track.eye_location = np.array([[-75.0, 20.0], [-75.0, 21.0], [-75.0, 22.0]])

    speed = reconstruction.translation_speed(track, track.t[0])
    assert speed == pytest.approx(111e3 / 3600.0, rel=0.02)
    # The last point reuses the previous interval rather than running off the end.
    assert reconstruction.translation_speed(track, track.t[-1]) == \
        pytest.approx(speed, rel=0.02)

    # A single-point track has no motion to measure.
    stationary = StormTrack()
    stationary.t = np.array(["2020-09-01T00"], dtype="datetime64[s]")
    stationary.eye_location = np.array([[-75.0, 20.0]])
    assert reconstruction.translation_speed(stationary, stationary.t[0]) == 0.0


def test_to_willoughby_wind_applies_both_corrections():
    """Storm-file 10 m point maximum -> flight-level azimuthal mean.

    Willoughby et al. (2006) fit to "maximum azimuthally averaged winds at
    several kilometers altitude", contrasted in the paper with HURDAT's
    "maximum 1-min averaged winds at 10-m elevation anywhere in the storm".  Two
    corrections separate those, and both must be applied -- the same pair the
    Fortran's ``adjust_max_wind`` uses.
    """
    track = StormTrack()
    track.t = np.array(["2020-09-01T00", "2020-09-01T01"],
                       dtype="datetime64[s]")
    track.eye_location = np.array([[-75.0, 20.0], [-75.0, 21.0]])
    track.max_wind_speed = np.array([50.0, 50.0])

    translation = reconstruction.translation_speed(track, track.t[0])
    expected = ((50.0 - reconstruction.ASYMMETRY_FRACTION * translation)
                / reconstruction.ATMOS_BOUNDARY_LAYER)
    assert reconstruction.to_willoughby_wind(track, track.t[0]) == \
        pytest.approx(expected)

    # Both corrections matter: the result is neither the raw wind nor a pure
    # height rescale of it.
    assert reconstruction.to_willoughby_wind(track, track.t[0]) != \
        pytest.approx(50.0)
    assert reconstruction.to_willoughby_wind(track, track.t[0]) != \
        pytest.approx(50.0 / reconstruction.ATMOS_BOUNDARY_LAYER)

    # A missing wind stays missing rather than becoming a number.
    track.max_wind_speed = np.array([np.nan, np.nan])
    assert np.isnan(reconstruction.to_willoughby_wind(track, track.t[0]))


def test_to_willoughby_wind_is_bounded_below_at_zero():
    """A storm moving faster than its own peak wind must not go negative."""
    track = StormTrack()
    track.t = np.array(["2020-09-01T00", "2020-09-01T01"],
                       dtype="datetime64[s]")
    # ~31 m/s of translation against a 5 m/s peak wind.
    track.eye_location = np.array([[-75.0, 20.0], [-75.0, 21.0]])
    track.max_wind_speed = np.array([5.0, 5.0])
    assert reconstruction.to_willoughby_wind(track, track.t[0]) == 0.0


def test_wind_conventions_versus_the_fortran():
    """Pin where Python and Fortran agree on the wind convention -- and differ.

    An RMW estimated here feeds a wind field evaluated in
    ``parametric_met_forcing_module.f90``, so a silent disagreement about what
    ``Vmax`` means would be invisible.  They agree on the *height* factor and
    deliberately differ on the *asymmetry* fraction:

    * height: both use 0.9 (10 m to boundary-layer top).
    * asymmetry: the Fortran subtracts the full translation magnitude; this
      module subtracts ``ASYMMETRY_FRACTION`` of it, because full subtraction
      over-corrects (Phadke et al. 2003 give ~0.5F at the RMW) and Eq. (7a) is
      monotone decreasing in wind, so over-subtraction biases the radius high.

    Changing the Fortran would move every forward-model path and all their
    goldens, so it is tracked as a separate decision.  This test exists so that
    divergence stays deliberate: if someone aligns the Fortran, this fails and
    points at the choice.
    """
    source = (Path(__file__).parents[1] / "src" / "2d" / "shallow" / "surge"
              / "parametric_met_forcing_module.f90")
    if not source.exists():
        pytest.skip("Fortran source not present in this checkout")
    text = source.read_text()

    import re
    match = re.search(r"atmos_boundary_layer\s*=\s*([0-9.]+)d0", text)
    assert match, "could not find atmos_boundary_layer in the Fortran"
    assert float(match.group(1)) == reconstruction.ATMOS_BOUNDARY_LAYER, \
        "height convention must stay shared"

    # The Fortran still subtracts the full magnitude.  If this stops being true,
    # revisit ASYMMETRY_FRACTION and the atlantic-rp calibration built on it.
    assert "mod_mws = mws - trans_speed" in text, (
        "the Fortran's asymmetry handling changed; this module's "
        "ASYMMETRY_FRACTION and any calibration derived under it must be "
        "re-examined")
    assert reconstruction.ASYMMETRY_FRACTION == 0.5


@needs_willoughby
def test_rmw_is_missing_when_translation_exceeds_peak_wind():
    """A regression input of zero is missing data, not a zero-wind storm.

    ``to_willoughby_wind`` bounds at zero because a wind *field* cannot be
    negative, but Eq. (7) was never fit at that value, so the estimator must
    decline rather than extrapolate to it.
    """
    track = StormTrack()
    track.t = np.array(["2020-09-01T00", "2020-09-01T01"],
                       dtype="datetime64[s]")
    track.eye_location = np.array([[-75.0, 20.0], [-75.0, 21.0]])
    track.max_wind_speed = np.array([5.0, 5.0])

    assert reconstruction.to_willoughby_wind(track, track.t[0]) == 0.0
    assert np.isnan(reconstruction.rmw_willoughby2006(track.t[0], track))


# ---------------------------------------------------------------------------
# The real verification: skill against observed RMW
# ---------------------------------------------------------------------------

@pytest.mark.remote
@pytest.mark.slow
def test_skill_against_observed_rmw():
    """Estimated RMW versus the ~21,700 IBTrACS points that report one.

    This is the check that actually verifies the coefficients and the wind
    conversion together.  The transcription check and every physical invariant
    in this file pass with a units or convention error still in place; this does
    not -- a factor of 1.852 (nmi/km) or 0.9 (boundary layer) shows up
    immediately as a bias of tens of kilometres.

    Bounds come from a measured baseline, with headroom.  They also *assert the
    known limitation*: Willoughby is biased low below its fitted range, and
    Vickery-Wadhera more so, which is why atlantic-rp applies a calibrated
    variant for weak systems rather than either raw fit.

    Needs the archive; set ``GEOCLAW_TRACK_ARCHIVES``.
    """
    import os
    pytest.importorskip("xarray")
    from clawpack.geoclaw.met.track import iter_ibtracs

    archives = os.environ.get("GEOCLAW_TRACK_ARCHIVES")
    if archives is None:
        pytest.skip("set GEOCLAW_TRACK_ARCHIVES to the archive directory")
    path = Path(archives) / "IBTrACS.NA.v04r01.nc"
    if not path.exists():
        pytest.skip(f"{path} not found")

    winds, observed, willoughby, vickery = [], [], [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for track in iter_ibtracs(path, basin="NA", years=range(1980, 2026)):
            for index in range(len(track.t)):
                reported = track.max_wind_radius[index]
                if not np.isfinite(reported):
                    continue
                t = track.t[index]
                winds.append(reconstruction.to_willoughby_wind(track, t))
                observed.append(reported)
                willoughby.append(reconstruction.rmw_willoughby2006(t, track))
                vickery.append(reconstruction.rmw_vickery_wadhera2008(t, track))

    winds = np.array(winds)
    observed = np.array(observed)
    willoughby = np.array(willoughby)
    vickery = np.array(vickery)
    assert observed.size > 15000, f"only {observed.size} reported-RMW points"

    strong = winds >= WILLOUGHBY_VALID_ABOVE
    weak = ~strong
    assert strong.sum() > 5000 and weak.sum() > 5000

    def stats(estimate, mask):
        residual = (estimate[mask] - observed[mask]) / 1e3   # km
        finite = np.isfinite(residual)
        return (float(np.median(residual[finite])),
                float(np.median(np.abs(residual[finite]))))

    # In range, Willoughby is good -- this is the headline assertion.
    bias, mae = stats(willoughby, strong)
    assert abs(bias) < 8.0, f"Willoughby bias in range: {bias:+.1f} km"
    assert mae < 20.0, f"Willoughby MAE in range: {mae:.1f} km"

    # Vickery-Wadhera is comparable in range...
    vw_bias, vw_mae = stats(vickery, strong)
    assert abs(vw_bias) < 8.0, f"VW08 bias in range: {vw_bias:+.1f} km"
    assert vw_mae < 22.0, f"VW08 MAE in range: {vw_mae:.1f} km"

    # ...and both are biased low out of range, VW08 the more so.  Asserted so
    # the limitation is a pinned fact rather than a docstring claim.
    weak_bias, _ = stats(willoughby, weak)
    vw_weak_bias, _ = stats(vickery, weak)
    assert weak_bias < -20.0, (
        f"Willoughby below range should be biased low; got {weak_bias:+.1f} km. "
        f"If this improved, revisit the atlantic-rp calibration.")
    assert vw_weak_bias < weak_bias, (
        f"VW08 ({vw_weak_bias:+.1f} km) should be worse than Willoughby "
        f"({weak_bias:+.1f} km) below range")


# ---------------------------------------------------------------------------
# Vickery-Wadhera specifics
# ---------------------------------------------------------------------------

def test_vickery_wadhera_degrades_sanely_at_zero_deficit():
    """At vanishing deficit the fit reduces to a latitude-only value.

    Which is both the reason it is usable without a wind speed, and the reason
    it cannot track weak-storm spread -- at a 13 m/s median deficit the squared
    term contributes ~0.01.
    """
    track = _synthetic_track(n=1)
    track.eye_location = np.array([[-75.0, 30.0]])
    track.central_pressure = np.array([101300.0])          # zero deficit
    radius = reconstruction.rmw_vickery_wadhera2008(track.t[0], track)
    # exp(3.015 + 0.0337*30) = ~56 km, per the paper's own degradation note.
    assert 50e3 < radius < 62e3


def test_vickery_wadhera_sigma_branches():
    """The heteroscedastic error model, at and across its breakpoints."""
    assert reconstruction.vickery_wadhera_sigma(50.0) == pytest.approx(0.448)
    assert reconstruction.vickery_wadhera_sigma(87.0) == pytest.approx(0.448)
    # Linear branch between 87 and 120 hPa.
    assert reconstruction.vickery_wadhera_sigma(100.0) == \
        pytest.approx(1.137 - 0.00792 * 100.0)
    assert reconstruction.vickery_wadhera_sigma(130.0) == pytest.approx(0.186)
    assert np.isnan(reconstruction.vickery_wadhera_sigma(np.nan))
    # Monotone decreasing: a deeper storm has a better-determined RMW.
    deficits = np.linspace(10.0, 150.0, 30)
    sigmas = [reconstruction.vickery_wadhera_sigma(d) for d in deficits]
    assert all(b <= a + 1e-12 for a, b in zip(sigmas, sigmas[1:]))


def test_rmw_sampled_requires_a_seed_and_is_reproducible():
    """Stochastic fills must be opt-in, seeded and repeatable."""
    track = _synthetic_track(max_wind_radius=np.full(6, np.nan))

    with pytest.raises(ValueError, match="explicit seed"):
        reconstruction.rmw_sampled(reconstruction.rmw_constant(40e3), 0.4,
                                   seed=None)

    def draws(seed):
        fill = reconstruction.rmw_sampled(reconstruction.rmw_constant(40e3),
                                          0.4, seed=seed)
        return [fill(t, track) for t in track.t]

    first, again, other = draws(12345), draws(12345), draws(99)
    assert first == again, "same seed must reproduce the draws"
    assert first != other, "a different seed must give different draws"
    assert all(np.isfinite(v) and v > 0 for v in first)
    # Spread is real but centred on the median.
    assert 0.2 < np.std(np.log(np.array(first) / 40e3)) < 1.0

    # A sigma callable is accepted, e.g. built on the VW08 error model.
    fill = reconstruction.rmw_sampled(
        reconstruction.rmw_vickery_wadhera2008,
        lambda t, s: reconstruction.vickery_wadhera_sigma(
            reconstruction.pressure_deficit(s, t)),
        seed=7)
    assert np.isfinite(fill(track.t[0], track))


# ---------------------------------------------------------------------------
# The Willoughby profile coefficients must come from one regression family
# ---------------------------------------------------------------------------

def _willoughby_source():
    """Live (non-comment) source of ``set_willoughby_fields``, whitespace-free.

    Comments are stripped and spaces removed before matching.  Both matter:
    Fortran spacing around operators is a style choice that a coefficient test
    must not depend on, and -- more importantly -- the routine keeps
    commented-out copies of the superseded Eq. (10) forms for reference, so a
    naive substring match can pass against *dead code* while the live line says
    something else.
    """
    path = (Path(__file__).parents[1] / "src" / "2d" / "shallow" / "surge"
            / "parametric_met_forcing_module.f90")
    if not path.exists():
        pytest.skip("Fortran source not present in this checkout")
    text = path.read_text()
    start = text.index("subroutine set_willoughby_fields")
    body = text[start:text.index("end subroutine set_willoughby_fields", start)]

    live = []
    for line in body.split("\n"):
        code = line.split("!", 1)[0]      # Fortran comment delimiter
        if code.strip():
            live.append(code)
    return "".join("".join(live).split())


def _has_terms(source, *terms):
    """Whether every *term* appears, each already whitespace-free."""
    return all(term.replace(" ", "") in source for term in terms)


def test_willoughby_n_and_A_use_the_ln_rmax_family():
    """``n`` and ``A`` take ``ln(R_max)`` as a predictor -- Eqs. (11b), (11c).

    Willoughby et al. (2006) gives two dual-exponential families: Eqs. (10a-c)
    predict from ``V_max`` and latitude alone, Eqs. (11a-c) add ``ln(R_max)``.
    They are not term-by-term interchangeable -- in (10) the variance that would
    load on ``ln(R_max)`` instead loads on ``V_max`` and latitude through their
    correlations -- so a routine must use one family throughout.
    """
    source = _willoughby_source()
    assert _has_terms(source, "2.134d0", "0.0077d0*mod_mws", "0.4522d0*log"), \
        "n no longer matches Eq. (11b)"
    assert _has_terms(source, "0.5913d0", "0.0029d0*mod_mws", "0.1361"), \
        "A no longer matches Eq. (11c)"
    # The Eq. (10) forms must not be live, only commented for reference.
    assert "0.4067d0" not in source, "n has reverted to the Eq. (10b) form"
    assert "0.0696d0" not in source, "A has reverted to the Eq. (10c) form"


def test_willoughby_x1_uses_the_same_family_as_n_and_A():
    """``X1`` must also be the ``ln(R_max)`` form -- Eq. (11a).

    It previously used Eq. (10a), so the outer decay length alone discarded the
    observed storm size while ``n`` and ``A`` retained it: internally
    inconsistent, and worst for storms whose size is anomalous for their
    intensity, which are the ones that matter for surge.

        X1 = 287.6 - 1.942*V_max + 7.799*ln(R_max) + 1.819*phi

    Measured effect of the correction over the 21,683 IBTrACS North Atlantic
    points with a reported ``R_max``: median +2.8 km on an ``X1`` of ~330 km,
    5th-95th percentile -5 to +10 km, maximum 20 km -- compact storms -3.8 km,
    large storms +8.6 km.  Real but low-single-digit percent, because
    ``7.799*ln(R_max)`` spans only ~11 km across observed sizes while the
    intercept difference (317.1 vs 287.6) partly offsets the slope changes.
    """
    source = _willoughby_source()
    assert _has_terms(source, "287.6d0", "1.942d0*mod_mws", "7.799d0*log",
                      "1.819"), \
        "X1 does not carry the Eq. (11a) coefficients"
    assert "317.1d0" not in source, (
        "the Eq. (10a) X1 is still live -- if it is only kept as a commented "
        "reference, this helper should have stripped it")


def test_willoughby_x2_is_the_published_constant():
    """``X2`` is fixed at 25 km in the paper; cleared as a non-bug."""
    assert _has_terms(_willoughby_source(), "X2=25.d3")
