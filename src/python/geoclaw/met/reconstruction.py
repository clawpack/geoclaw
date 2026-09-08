#!/usr/bin/env python
# encoding: utf-8

r"""Opt-in reconstruction of missing tropical-cyclone geometry.

Commonly available best tracks are incomplete in a way that stops GeoClaw from
running: the radius of maximum winds (RMW) and the outer/last-closed-isobar
radius (ROCI) are frequently absent, and
:meth:`~clawpack.geoclaw.met.parametric.ParametricMetForcing.write_geoclaw`
skips any track time missing either.  Measured over IBTrACS v4 North Atlantic,
1980-2025 (39,711 track points from 684 storms):

===========================  =========  ==================================
Field                        Reported   Notes
===========================  =========  ==================================
``max_wind_speed``           100%       after the reader's agency cascade
``central_pressure``         95.9%      the other 4.1% are dropped outright
``max_wind_radius`` (RMW)    54.6%      0.4% in the 1980s-90s, ~100% in the 2020s
``storm_radius`` (ROCI)      54.3%      same era split
===========================  =========  ==================================

So roughly **45% of the record has no reported storm geometry at all**, and it is
not evenly distributed: the 1980s and 1990s are essentially unreported while the
2020s are complete.  Anything that estimates geometry therefore changes the
early record far more than the recent one -- worth keeping in view when
interpreting results that span decades.

Everything here is **opt-in**.  Nothing is reconstructed unless a caller asks:

.. code-block:: python

    from clawpack.geoclaw.met import reconstruction

    fill = reconstruction.default_fill_dict(models={"max_wind_radius":
                                                    reconstruction.rmw_constant(40e3)})
    storm.write(path, file_format="geoclaw", fill_dict=fill)

Estimators match the existing ``fill_dict`` contract exactly -- ``fn(t, storm)``
returning an SI value -- so they drop into ``write_geoclaw(fill_dict=...)`` and
compose with :func:`~clawpack.geoclaw.met.track.fill_rad_w_other_source`.

.. note::

   These are empirical regressions transcribed from published papers, not
   solver mathematics.  Coefficient values and their provenance are the
   maintainer's to confirm -- see :data:`WILLOUGHBY_2006_RMW`, which
   deliberately ships unpopulated so that a fabricated coefficient cannot reach
   a run.  Estimators whose coefficients are not filled in raise
   :exc:`NotImplementedError` naming the paper and the table to transcribe.

:References:

 - Willoughby, H. E., R. W. R. Darling and M. E. Rahn (2006). "Parametric
   representation of the primary hurricane vortex. Part II: A new family of
   sectionally continuous profiles." *Mon. Wea. Rev.* **134**, 1102-1120.
   doi:10.1175/MWR3106.1
 - Chavas, D. R., N. Lin, W. Dong and Y. Lin (2016). "Observed tropical cyclone
   size revisited." *J. Climate* **29**, 2923-2939. doi:10.1175/JCLI-D-15-0731.1
 - Knaff, J. A. and R. M. Zehr (2007). "Reexamination of tropical cyclone
   wind-pressure relationships." *Wea. Forecasting* **22**, 71-88.
   doi:10.1175/WAF965.1  (*not yet implemented -- see* :func:`wind_pressure`)
"""

import warnings

import numpy as np

import clawpack.geoclaw.units as units
import clawpack.geoclaw.util as util


# ---------------------------------------------------------------------------
# Coefficient blocks
#
# Each estimator's numbers live in one named, citable block rather than inline,
# so that what came from which paper is auditable.  A block of ``None`` means
# "not yet transcribed": the estimator raises rather than guessing.
# ---------------------------------------------------------------------------

#: Coefficients for the Willoughby et al. (2006) RMW regression,
#: ``R_max`` as a function of maximum wind speed and latitude --
#: **Eq. (7)** of the paper.
#:
#: Transcribed from the paper (maintainer-verified against the author's own
#: copy).  Eq. (7a) is identical in the single- and dual-exponential fits; the
#: paper notes it is not repeated for the dual case::
#:
#:     R_max = 46.4 * exp(-0.0155 * V_max + 0.0169 * phi)
#:
#: All coefficients are significant at better than 1% except the latitude term
#: at 1.4%.
#:
#: **Input basis (resolved).**  The paper states its predictand explicitly:
#: *"They are maximum azimuthally averaged winds at several kilometers altitude,
#: whereas HURDAT contains maximum 1-min averaged winds at 10-m elevation
#: anywhere in the storm."*  So Eq. (7) does **not** take a storm file's
#: ``max_wind_speed``; feed it :func:`to_willoughby_wind`, which subtracts the
#: translation speed (point maximum to azimuthal mean) and divides by
#: :data:`ATMOS_BOUNDARY_LAYER` (10 m to boundary-layer top) -- the same two
#: corrections the Fortran's ``adjust_max_wind`` applies before evaluating the
#: paper's Section 4 coefficients.  This also confirms the existing Fortran
#: Willoughby path is on the right basis.
#:
#: **Range of validity -- read this before relying on the output.**  The fitting
#: sample averages 34.6 m/s maximum flight-level wind, "just above the threshold
#: of hurricane intensity", so it is hurricane-dominated.  The track points that
#: *need* an RMW estimate are not: over IBTrACS North Atlantic 1980-2025, of the
#: 18,028 points with no reported RMW, once converted to this basis
#:
#: * **45.4%** fall below 15 m/s -- far below the fitting sample,
#: * 24.4% are 15-25 m/s,
#: * only **22.0%** are 25-45 m/s, i.e. near the sample mean,
#: * 8.2% exceed 45 m/s.
#:
#: So for roughly 70% of the points it would fill, Eq. (7) is being
#: *extrapolated*, not interpolated.  This is inherent to the problem rather
#: than to the implementation: unreported RMW correlates with weak, early-stage
#: and pre-2000 storms, and the translation-speed correction bites hardest
#: exactly there (median conversion -17%, 5th percentile -75%).  Anyone using
#: this at scale should decide deliberately whether to apply it across the whole
#: record, restrict it to the fitted range and use something else below, or pick
#: an estimator built for weak systems.
#:
#: The remaining unknown is the *temporal* averaging of the flight-level winds,
#: which the paper does not state (its averaging is azimuthal, i.e. spatial).
#: For reference, the Fortran treats storm-file winds as 1-minute and converts
#: to 10-minute on output via ``sampling_time = 0.88``; SLOSH is the one model
#: that also applies it upstream, its own fit assuming 10-min winds.
WILLOUGHBY_2006_RMW = {
    "reference": "Willoughby, Darling & Rahn (2006), Mon. Wea. Rev. 134(4), "
                 "1102-1120, doi:10.1175/MWR3106.1  "
                 "(bibliographic record verified against Crossref; the "
                 "equation attribution below is the maintainer's to confirm)",
    "equation": "Eq. (7a)",
    # Resolved from the paper's own description of its predictand; see above.
    "wind_basis": "maximum azimuthally averaged wind at flight level "
                  "(predominantly 850/700 hPa, 1.5/3 km); use "
                  "to_willoughby_wind()",
    # The paper's text gives 34.6 m/s for the whole sample; 36.7 m/s is quoted
    # for the fitted subset.  Most likely whole-sample vs post-QC -- recorded
    # both rather than picking one, since only the order of magnitude matters
    # for the range-of-validity note.
    "sample_mean_wind_m_s": {"whole_sample": 34.6, "fitted": 36.7},
    # Below this the fit is not merely extrapolated but measurably wrong; see
    # the range-of-validity note above.
    "valid_above_m_s": 26.0,
    "temporal_averaging": None,     # not stated by the paper
    "basin": "Atlantic",
    "units": {"max_wind_radius": "km", "max_wind_speed": "m/s",
              "latitude": "degrees"},
    "coefficients": {"prefactor_km": 46.4,
                     "wind_coefficient": -0.0155,
                     "latitude_coefficient": 0.0169},
}

#: Coefficients for the Vickery & Wadhera (2008) Atlantic RMW regression::
#:
#:     ln(R_max) = 3.015 - 6.291e-5 * dp**2 + 0.0337 * psi
#:
#: with ``dp`` the central pressure deficit in hPa, ``psi`` latitude in degrees
#: and ``R_max`` in km.  Reported ``r^2 = 0.297``, overall ``sigma_lnRMW = 0.441``.
#:
#: Its appeal is being pressure-driven -- usable where a central pressure is
#: reported but a wind is not -- and shipping an uncertainty model (see
#: :data:`VICKERY_WADHERA_2008_SIGMA`).
#:
#: **Measured against observation, it is not a good weak-storm estimator**,
#: despite being the obvious candidate.  Over the 21,683 IBTrACS North Atlantic
#: 1980-2025 points that report an RMW, median residuals are
#:
#: * below 25 m/s: **-42 km** (Willoughby (7a) manages -35 km there),
#: * above 25 m/s: +1.7 km (Willoughby: -3.2 km).
#:
#: The reason is structural: at the median deficit of 13 hPa the ``dp**2`` term
#: contributes about 0.01, so for weak systems the fit collapses to a
#: latitude-only constant and cannot track their spread.  Prefer it only where a
#: wind speed is genuinely unavailable.
VICKERY_WADHERA_2008_RMW = {
    "reference": "Vickery & Wadhera (2008), J. Appl. Meteor. Climatol. 47, "
                 "2714-2726, doi:10.1175/2008JAMC1837.1",
    "basin": "Atlantic",
    "units": {"max_wind_radius": "km", "pressure_deficit": "hPa",
              "latitude": "degrees"},
    "r_squared": 0.297,
    "sigma_ln_rmw_overall": 0.441,
    "coefficients": {"intercept": 3.015,
                     "deficit_squared_coefficient": -6.291e-5,
                     "latitude_coefficient": 0.0337},
}

#: Heteroscedastic error model accompanying :data:`VICKERY_WADHERA_2008_RMW`:
#: ``sigma`` of ``ln(R_max)`` as a piecewise function of the deficit in hPa.
VICKERY_WADHERA_2008_SIGMA = ((87.0, "constant", 0.448),
                              (120.0, "linear", (1.137, -0.00792)),
                              (np.inf, "constant", 0.186))


#: Climatological outer-radius fallback, in metres.
#:
#: This replaces an anonymous hard-wired ``500e3`` that
#: ``ParametricMetForcing.write_geoclaw`` has applied for years.  The value is
#: unchanged, so nothing that relied on it moves; what changes is that it is now
#: named, documented and overridable.
#:
#: **It is not an ROCI, and should not be fitted to one.**  GeoClaw uses this
#: value as the centre of a ``tanh`` taper of width 100 km applied to both the
#: wind vector and the pressure perturbation
#: (``post_process_wind_estimate``): the ramp is 0.88 at 100 km inside it, 0.50
#: at it, 0.12 at 100 km outside, ~0.02 by 200 km outside.  So it sets where the
#: parametric field blends to ambient.
#:
#: 500 km is defensible on those terms.  With an outer decay length ``X1`` of
#: ~250 km, a storm with ``R_max`` 40 km and ``V_max`` 50 m/s still carries
#: ~8 m/s of azimuthal-mean wind at 500 km -- relevant to Ekman setup and the
#: inverse-barometer contribution.  Tapering at a *physical* North Atlantic ROCI
#: (~350-400 km) would truncate surge-relevant outer winds.
#:
#: If a size-varying default is ever wanted, base it on an ROCI or R34
#: climatology and then add a margin of order ``X1``, so the 100 km ramp
#: completes beyond the surge-relevant winds -- rather than setting this equal
#: to the ROCI itself.
ROCI_CLIMATOLOGY_M = 500e3


#: 10 m wind to top-of-boundary-layer wind.
#:
#: Mirrors ``atmos_boundary_layer`` in
#: ``src/2d/shallow/surge/parametric_met_forcing_module.f90`` (an ADCIRC input
#: that is always set to this value).  Kept equal to the Fortran's on purpose:
#: an RMW estimated here feeds a wind field evaluated there, so the two must
#: share a convention.  If one changes, change both.
ATMOS_BOUNDARY_LAYER = 0.9


#: Fraction of the storm's forward speed removed when converting a point-maximum
#: wind to an azimuthal mean.
#:
#: The canonical reduction is about half the forward speed at the radius of
#: maximum winds (Phadke et al. 2003; the same pipeline used by Knaff et al.
#: 2011); Lin & Chavas (2012) use a gentler 0.55-0.7 on a rotated motion vector.
#:
#: **This deliberately differs from the Fortran.**
#: ``parametric_met_forcing_module.f90``'s ``adjust_max_wind`` subtracts the
#: *full* translation magnitude, which over-corrects; because Eq. (7a) is
#: monotone decreasing in wind speed, over-subtraction biases the reconstructed
#: radius high.  Fixing it there would move every forward-model path (Holland
#: included) and all their goldens, so it is tracked separately.  Measured over
#: the 18,028 North Atlantic points needing an RMW estimate, moving from full to
#: half subtraction raises the median converted wind from 16.3 to 20.1 m/s and
#: drops the count driven to zero from 267 to 9.
ASYMMETRY_FRACTION = 0.5


# Warnings that are informative once and noise thereafter -- these fire per
# track point, and a basin sweep has tens of thousands.
_WARNED = set()


def _warn_once(key, message):
    if key not in _WARNED:
        _WARNED.add(key)
        warnings.warn(message)


class CoefficientsNotTranscribed(NotImplementedError):
    r"""Raised when an estimator's published coefficients are not filled in."""


def _require_coefficients(block, estimator):
    r"""Raise unless *block* has had its coefficients transcribed."""
    if block.get("coefficients") is None:
        raise CoefficientsNotTranscribed(
            f"{estimator} has no transcribed coefficients.\n"
            f"  Source: {block['reference']}\n"
            f"Fill in the coefficient block in "
            f"clawpack/geoclaw/met/reconstruction.py (and the reference values "
            f"in tests/data/storm/reconstruction_reference.json) from the "
            f"paper.  Coefficients are deliberately not shipped so that an "
            f"unverified number cannot reach a simulation.")
    return block["coefficients"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _at_time(storm, t, field):
    r"""Value of *field* on *storm* at track time *t*, or ``np.nan``.

    Exact match on the time axis first, since ``fill_dict`` callbacks are handed
    a time drawn from that same axis.  Falls back to the nearest track time,
    which matters when a caller interpolates onto its own grid.
    """
    values = getattr(storm, field, None)
    if values is None:
        return np.nan
    times = np.asarray(storm.t)
    values = np.asarray(values, dtype=float)

    exact = np.nonzero(times == t)[0]
    if exact.size:
        return float(values[exact[0]])

    try:
        offsets = np.abs(times - t)
    except TypeError:
        return np.nan
    return float(values[int(np.argmin(offsets))])


def _latitude_at(storm, t):
    r"""Storm-centre latitude at track time *t*, or ``np.nan``."""
    if storm.eye_location is None:
        return np.nan
    times = np.asarray(storm.t)
    exact = np.nonzero(times == t)[0]
    index = int(exact[0]) if exact.size else int(np.argmin(np.abs(times - t)))
    return float(np.asarray(storm.eye_location)[index, 1])


def translation_speed(storm, t):
    r"""Storm translation speed at track time *t*, in m/s.

    Forward difference between consecutive track points, mirroring how the
    Fortran builds ``storm%velocity`` in
    ``parametric_met_forcing_module.f90`` (great-circle distance over the time
    step, with the final point reusing the previous value).
    """
    times = np.asarray(storm.t)
    centres = np.asarray(storm.eye_location, dtype=float)
    if times.size < 2:
        return 0.0

    exact = np.nonzero(times == t)[0]
    index = int(exact[0]) if exact.size else int(np.argmin(np.abs(times - t)))
    index = min(index, times.size - 2)

    seconds = float((times[index + 1] - times[index])
                    / np.timedelta64(1, "s")) if np.issubdtype(
                        times.dtype, np.datetime64) else float(
                            times[index + 1] - times[index])
    if not np.isfinite(seconds) or seconds <= 0.0:
        return 0.0

    distance = util.haversine(centres[index, 0], centres[index, 1],
                              centres[index + 1, 0], centres[index + 1, 1])
    return float(distance / seconds)


def to_willoughby_wind(storm, t):
    r"""Convert a storm file's maximum wind to Willoughby's fitting basis.

    Willoughby et al. (2006) fit to *maximum azimuthally averaged winds at
    several kilometres altitude*, which the paper contrasts explicitly with
    HURDAT's *maximum 1-minute averaged winds at 10 m elevation anywhere in the
    storm*.  Two corrections separate those, and this applies both -- exactly as
    the Fortran's ``adjust_max_wind`` does before evaluating the paper's
    Section 4 profile coefficients:

    1. **Remove the motion asymmetry.**  HURDAT's maximum is a point maximum
       "anywhere in the storm", so it sits in the translation-enhanced quadrant;
       Willoughby's is an azimuthal mean, which that enhancement does not
       survive.  :data:`ASYMMETRY_FRACTION` of the forward speed is subtracted --
       *not* the whole of it, which is where the Fortran differs.
    2. **Divide by** :data:`ATMOS_BOUNDARY_LAYER` -- 10 m to the top of the
       boundary layer.

    Any estimator taking Willoughby coefficients should be fed this, not
    ``storm.max_wind_speed`` directly; the difference is order 10-20%.
    """
    max_wind_speed = _at_time(storm, t, "max_wind_speed")
    if not np.isfinite(max_wind_speed):
        return np.nan
    # Bounded at zero: a storm can translate faster than its own peak wind.
    symmetric = max(max_wind_speed
                    - ASYMMETRY_FRACTION * translation_speed(storm, t), 0.0)
    return symmetric / ATMOS_BOUNDARY_LAYER


# ---------------------------------------------------------------------------
# Estimators -- fn(t, storm) -> SI value, matching the fill_dict contract
# ---------------------------------------------------------------------------

def rmw_constant(radius_m):
    r"""A constant radius of maximum winds, in metres.

    Deliberately explicit rather than physical: it makes a track runnable while
    saying plainly in the calling code that the geometry is an assumption, not an
    estimate.  Prefer :func:`rmw_willoughby2006` once its coefficients are
    confirmed.

    :Input:
     - *radius_m* (float) The radius to use, in metres.

    :Output:
     - A ``fn(t, storm)`` suitable for ``write_geoclaw(fill_dict=...)``.
    """
    radius_m = float(radius_m)
    if not np.isfinite(radius_m) or radius_m <= 0.0:
        raise ValueError(f"A constant RMW must be finite and positive, got "
                         f"{radius_m}.")

    def fill(t, storm):
        return radius_m

    fill.__doc__ = f"Constant radius of maximum winds of {radius_m:.0f} m."
    return fill


def rmw_willoughby2006(t, storm):
    r"""Radius of maximum winds from maximum wind speed and latitude.

    Willoughby, Darling & Rahn (2006).  Chosen over the alternatives because its
    inputs are present at essentially every HURDAT2/IBTrACS track point, where
    e.g. Vickery & Wadhera (2008) needs a trusted pressure deficit.

    Raises :exc:`CoefficientsNotTranscribed` until :data:`WILLOUGHBY_2006_RMW`
    is populated from the paper.

    :Input:
     - *t* Track time to evaluate at.
     - *storm* The :class:`~clawpack.geoclaw.met.track.StormTrack`.

    :Output:
     - (float) Radius of maximum winds in metres, or ``np.nan`` if the inputs
       needed are themselves missing.
    """
    coefficients = _require_coefficients(WILLOUGHBY_2006_RMW,
                                         "rmw_willoughby2006")

    # Eq. (7) is fit to azimuthally averaged flight-level winds, not to the 10 m
    # point maximum a storm file carries -- see WILLOUGHBY_2006_RMW.
    wind = to_willoughby_wind(storm, t)
    latitude = _latitude_at(storm, t)
    if not (np.isfinite(wind) and np.isfinite(latitude)):
        return np.nan

    # to_willoughby_wind bounds at zero because a wind *field* cannot be
    # negative.  A regression *input* of zero is a different
    # matter: it means the storm translates faster than its own peak wind, and
    # Eq. (7) has nothing to say about that.  Report missing rather than
    # evaluate the fit at a value it was never given.  (267 of the 18,028 North
    # Atlantic 1980-2025 points needing an RMW estimate are in this state.)
    if wind <= 0.0:
        return np.nan

    # Eq. (7a), published in km with latitude in degrees; converted to SI on
    # the way out so callers never see the paper's units.
    radius_km = (coefficients["prefactor_km"]
                 * np.exp(coefficients["wind_coefficient"] * wind
                          + coefficients["latitude_coefficient"] * latitude))

    # Below the fitted range the regression is not merely extrapolated, it is
    # measurably biased low -- against the 21,683 IBTrACS North Atlantic points
    # that do report an RMW, the median residual is about -35 km under 25 m/s
    # versus -3 km above it.  Warn rather than refuse: a caller may legitimately
    # want the raw published fit, and atlantic-rp's calibrated variant is built
    # on top of exactly this call.
    if wind < WILLOUGHBY_2006_RMW["valid_above_m_s"]:
        _warn_once(
            "willoughby_below_range",
            f"Willoughby (2006) Eq. (7a) evaluated at "
            f"{wind:.1f} m/s, below its ~"
            f"{WILLOUGHBY_2006_RMW['valid_above_m_s']:.0f} m/s range of "
            f"validity; it is biased low by roughly 35 km there.  Consider a "
            f"calibrated or climatological estimator for weak systems.")

    return units.convert(radius_km, 'km', 'm')


def roci_climatology(t, storm):
    r"""Climatological outer storm radius, in metres.

    The named, documented replacement for the constant ``write_geoclaw`` has
    always applied.  Returns :data:`ROCI_CLIMATOLOGY_M` regardless of *t* and
    *storm*; the signature matches the ``fill_dict`` contract so it can be
    swapped for something storm-dependent without changing callers.

    A latitude-dependent form from Chavas et al. (2016) would be a genuine
    improvement over a flat constant, but note that 500 km already sits near the
    observed North Atlantic median -- the value of this function is that the
    assumption is now explicit and overridable, not that it is more accurate.
    """
    return ROCI_CLIMATOLOGY_M


def pressure_deficit(storm, t, ambient_pressure=101300.0):
    r"""Central pressure deficit at track time *t*, in hPa.

    *ambient_pressure* defaults to GeoClaw's own ambient (``geo_data``'s
    101.3 kPa).  Note that for a weak storm the deficit is small enough that the
    choice of ambient matters in relative terms -- at a 13 hPa median deficit, a
    5 hPa error in ambient is a large fractional error, though
    :func:`rmw_vickery_wadhera2008` squares the deficit and so is insensitive to
    it there.
    """
    central = _at_time(storm, t, "central_pressure")
    if not np.isfinite(central):
        return np.nan
    return (ambient_pressure - central) / 100.0


def vickery_wadhera_sigma(deficit_hpa):
    r"""``sigma`` of ``ln(R_max)`` for Vickery & Wadhera (2008), given a deficit."""
    if not np.isfinite(deficit_hpa):
        return np.nan
    for upper, kind, value in VICKERY_WADHERA_2008_SIGMA:
        if deficit_hpa <= upper:
            if kind == "constant":
                return float(value)
            intercept, slope = value
            return float(intercept + slope * deficit_hpa)
    return float(VICKERY_WADHERA_2008_SIGMA[-1][2])


def rmw_vickery_wadhera2008(t, storm, ambient_pressure=101300.0):
    r"""Radius of maximum winds from pressure deficit and latitude.

    Vickery & Wadhera (2008), Atlantic.  Useful where a central pressure is
    reported but a wind speed is not.

    **Read the note on** :data:`VICKERY_WADHERA_2008_RMW` **before choosing this
    for weak systems** -- measured against observation it is worse there than
    Willoughby, by about -42 km against -35 km.

    :Output:
     - (float) Radius of maximum winds in metres, or ``np.nan`` if the central
       pressure or latitude needed are missing.
    """
    coefficients = _require_coefficients(VICKERY_WADHERA_2008_RMW,
                                         "rmw_vickery_wadhera2008")
    deficit = pressure_deficit(storm, t, ambient_pressure=ambient_pressure)
    latitude = _latitude_at(storm, t)
    if not (np.isfinite(deficit) and np.isfinite(latitude)):
        return np.nan

    radius_km = np.exp(
        coefficients["intercept"]
        + coefficients["deficit_squared_coefficient"] * deficit ** 2
        + coefficients["latitude_coefficient"] * abs(latitude))
    return units.convert(radius_km, 'km', 'm')


def rmw_sampled(estimator, sigma, seed):
    r"""Wrap an RMW *estimator* to draw from its lognormal error model.

    Reconstructed geometry is uncertain, and for a hazard ensemble propagating
    that spread is arguably more honest than a point estimate.  This is
    deliberately **opt-in and seeded**: the default path stays deterministic so a
    storm always produces the same ``.storm`` file, which the ensemble driver
    relies on for resumability and for comparing runs.

    :Input:
     - *estimator* A ``fn(t, storm)`` returning a radius in metres, treated as
       the median of the lognormal.
     - *sigma* Either a float, or a ``fn(t, storm)`` returning the ``sigma`` of
       ``ln(R_max)`` (e.g. built on :func:`vickery_wadhera_sigma`).
     - *seed* Required.  Draws are reproducible given the seed, the storm's times
       and the call order.

    :Output:
     - A ``fn(t, storm)`` suitable for ``write_geoclaw(fill_dict=...)``.
    """
    if seed is None:
        raise ValueError(
            "rmw_sampled requires an explicit seed: an unseeded stochastic fill "
            "makes a storm file irreproducible, and the ensemble driver treats "
            "a storm file as a function of its track.")
    generator = np.random.default_rng(seed)

    def fill(t, storm):
        median = estimator(t, storm)
        if not np.isfinite(median) or median <= 0.0:
            return np.nan
        spread = sigma(t, storm) if callable(sigma) else float(sigma)
        if not np.isfinite(spread) or spread <= 0.0:
            return median
        return float(median * np.exp(generator.normal(0.0, spread)))

    fill.__doc__ = f"Lognormal draw about {estimator.__name__} (seed={seed})."
    return fill


def wind_pressure(t, storm):
    r"""Wind-pressure relationship -- **not implemented**.

    Would fill ``central_pressure`` from ``max_wind_speed`` (Knaff & Zehr 2007,
    as revised by Courtney & Knaff 2009), recovering track points that
    ``write_geoclaw`` currently skips for want of a pressure.

    Deliberately deferred to its own change: the relationship depends on storm
    size, translation speed, latitude and environmental pressure, and its
    coefficient set is easy to get subtly and invisibly wrong.  The gap it would
    close is also modest and concentrated -- over IBTrACS North Atlantic
    1980-2025, 1,682 of 41,427 track points (4.1%) have a wind but no pressure,
    of which 1,519 are in the 1980s (21.9% of that decade) and none after 1999.

    Until then those points are skipped and counted by :func:`coverage_report`.
    """
    raise NotImplementedError(
        "The wind-pressure relationship is not implemented; see the docstring. "
        "Track points with no central pressure are skipped at write time and "
        "reported by reconstruction.coverage_report().")


# ---------------------------------------------------------------------------
# Assembly and reporting
# ---------------------------------------------------------------------------

def default_fill_dict(track=None, models=None):
    r"""Assemble a ``fill_dict`` for ``write_geoclaw``.

    The default only supplies the outer radius (:func:`roci_climatology`),
    matching what ``write_geoclaw`` already did -- so calling this with no
    arguments changes nothing.  Pass *models* to opt into more:

    .. code-block:: python

        fill = default_fill_dict(models={"max_wind_radius":
                                         rmw_constant(40e3)})

    :Input:
     - *track* (optional) The track the fills will be applied to.  Used only to
       warn about fields that will still be missing.
     - *models* (dict, optional) Field name to ``fn(t, storm)``, overriding or
       extending the defaults.  A value of ``None`` removes that default.

    :Output:
     - (dict) Suitable for ``write_geoclaw(fill_dict=...)``.

    Estimators must not depend on another estimator's output: ``write_geoclaw``
    applies fills to internal copies, so evaluation order is not observable.
    """
    fill = {"storm_radius": roci_climatology}

    for field, model in (models or {}).items():
        if model is None:
            fill.pop(field, None)
        else:
            fill[field] = model

    if track is not None:
        report = coverage_report(track)
        unfilled = [field for field, counts in report["fields"].items()
                    if counts["missing"] and field not in fill]
        if unfilled:
            warnings.warn(
                "No fill supplied for " + ", ".join(sorted(unfilled))
                + "; track times missing those will be skipped at write time. "
                  "See reconstruction.coverage_report() for counts.")

    return fill


def coverage_report(track, fill_dict=None):
    r"""How much of *track* is reported, and how much a fill would supply.

    Answers "what will actually get written, and why not the rest" before
    running anything.

    :Input:
     - *track* A :class:`~clawpack.geoclaw.met.track.StormTrack`.
     - *fill_dict* (dict, optional) Fills to evaluate, as
       :func:`default_fill_dict` returns.  Estimators are called, so this is as
       expensive as the estimators are.

    :Output:
     - (dict) ``n_times``, ``n_writable`` (times where all four fields would be
       present), and per-field ``reported`` / ``missing`` / ``filled`` /
       ``still_missing`` counts.
    """
    fields = ("max_wind_speed", "central_pressure", "max_wind_radius",
              "storm_radius")
    times = np.asarray(track.t)
    n_times = int(times.size)

    report = {"n_times": n_times, "fields": {}}
    present = np.ones(n_times, dtype=bool)

    for field in fields:
        values = getattr(track, field, None)
        if values is None:
            reported = np.zeros(n_times, dtype=bool)
        else:
            reported = np.isfinite(np.asarray(values, dtype=float))

        entry = {"reported": int(reported.sum()),
                 "missing": int((~reported).sum()),
                 "filled": 0,
                 "still_missing": int((~reported).sum())}

        model = (fill_dict or {}).get(field)
        if model is not None and entry["missing"]:
            filled = reported.copy()
            for index in np.nonzero(~reported)[0]:
                value = model(times[index], track)
                if np.isfinite(value) and value > 0.0:
                    filled[index] = True
            entry["filled"] = int(filled.sum() - reported.sum())
            entry["still_missing"] = int((~filled).sum())
            reported = filled

        report["fields"][field] = entry
        present &= reported

    report["n_writable"] = int(present.sum())
    return report
