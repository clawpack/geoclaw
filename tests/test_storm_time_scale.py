"""
Self-contained test for the storm_time_scale time-axis arithmetic.

No GeoClaw runtime or test framework is required — only numpy and assert.
Mirrors the Fortran arithmetic in data_storm_module.f90 case(2):

    t_ref = t[0]
    t_scaled[i] = t_ref + round(scale * (t[i] - t_ref))
"""

import numpy as np


def apply_time_scale(t_raw, scale):
    """Apply storm_time_scale to a raw time array.

    Parameters
    ----------
    t_raw : array-like of float
        Times in any consistent unit (e.g. seconds since epoch).
    scale : float
        storm_time_scale value.

    Returns
    -------
    numpy.ndarray of int
        Scaled integer times matching Fortran nint() rounding.
    """
    t = np.asarray(t_raw, dtype=np.float64)
    t_ref = t[0]
    return np.rint(t_ref + scale * (t - t_ref)).astype(np.int64)


def _hourly_times(n=10):
    """Return n hourly timestamps as float64 seconds since Unix epoch."""
    epoch = np.datetime64("1970-01-01T00:00:00", "s")
    start = np.datetime64("2012-08-28T00:00:00", "s")
    t0 = float((start - epoch) / np.timedelta64(1, "s"))
    return np.array([t0 + i * 3600.0 for i in range(n)])


def test_scale_one_is_noop():
    t = _hourly_times(10)
    scaled = apply_time_scale(t, 1.0)
    expected = np.rint(t).astype(np.int64)
    assert np.array_equal(scaled, expected), f"scale=1.0 changed values: {scaled - expected}"


def test_scale_two_doubles_duration():
    t = _hourly_times(10)
    scaled = apply_time_scale(t, 2.0)
    # Duration from first to last should be doubled
    orig_duration = t[-1] - t[0]
    scaled_duration = scaled[-1] - scaled[0]
    assert abs(scaled_duration - 2.0 * orig_duration) < 1, (
        f"Expected duration {2*orig_duration}, got {scaled_duration}"
    )


def test_anchor_invariant():
    t = _hourly_times(10)
    t_ref = np.rint(t[0]).astype(np.int64)
    for scale in [0.5, 1.0, 2.0, 3.7]:
        scaled = apply_time_scale(t, scale)
        assert scaled[0] == t_ref, (
            f"scale={scale}: anchor changed from {t_ref} to {scaled[0]}"
        )


def test_scale_half_compresses_duration():
    t = _hourly_times(10)
    scaled = apply_time_scale(t, 0.5)
    orig_duration = t[-1] - t[0]
    scaled_duration = scaled[-1] - scaled[0]
    assert abs(scaled_duration - 0.5 * orig_duration) < 1, (
        f"Expected duration {0.5*orig_duration}, got {scaled_duration}"
    )


# ---------------------------------------------------------------------------
# Model-storm time scale (Cap. 4A)
# Mirrors set_storm_fields in storm_module.f90:
#     t_eff = landfall + (t - landfall) / storm_time_scale
# Division (not multiplication) makes scale>1 = slower, matching the data
# storm: a given track state is reached at a later sim time.
# ---------------------------------------------------------------------------

def model_time_eff(t, landfall, scale):
    return landfall + (t - landfall) / scale


def test_model_scale_one_is_noop():
    for t in (-3600.0, 0.0, 7200.0):
        assert model_time_eff(t, 0.0, 1.0) == t


def test_model_scale_anchors_at_landfall():
    landfall = 1000.0
    for scale in (0.5, 2.0, 3.7):
        assert model_time_eff(landfall, landfall, scale) == landfall


def test_model_scale_direction_matches_data_storm():
    # scale>1 = slower: at a fixed sim time past landfall, the effective
    # track time is closer to landfall (storm has progressed less) than the
    # unscaled case.
    landfall, t = 0.0, 7200.0
    t_eff_slow = model_time_eff(t, landfall, 2.0)   # slower
    t_eff_norm = model_time_eff(t, landfall, 1.0)
    assert t_eff_slow < t_eff_norm
    # The data-storm relabeling (stretch by scale) and the model-storm lookup
    # (divide by scale) are inverse maps: composing them is the identity.
    data_sim_time = landfall + 2.0 * (t_eff_slow - landfall)
    assert abs(data_sim_time - t) < 1e-9


# ---------------------------------------------------------------------------
# Temporal onset/cutoff ramp (Cap. 4B)
# Mirrors temporal_ramp() in storm_module.f90 (raised cosine, 0 at the
# window edge -> 1 inside over the ramp width; 0 width disables that side).
# ---------------------------------------------------------------------------

def temporal_ramp(t, t0, tfinal, t_ramp_on, t_ramp_off):
    onset = 1.0
    if t_ramp_on > 0.0 and t < t0 + t_ramp_on:
        onset = 0.5 * (1.0 - np.cos(np.pi * max(0.0, t - t0) / t_ramp_on))
    cutoff = 1.0
    if t_ramp_off > 0.0 and t > tfinal - t_ramp_off:
        cutoff = 0.5 * (1.0 - np.cos(np.pi * max(0.0, tfinal - t) / t_ramp_off))
    return onset * cutoff


def test_temporal_ramp_disabled_is_unity():
    for t in (0.0, 500.0, 1000.0):
        assert temporal_ramp(t, 0.0, 1000.0, 0.0, 0.0) == 1.0


def test_temporal_ramp_onset():
    t0, tfinal, w = 0.0, 1000.0, 100.0
    assert abs(temporal_ramp(t0, t0, tfinal, w, 0.0) - 0.0) < 1e-12      # off at t0
    assert abs(temporal_ramp(t0 + w / 2, t0, tfinal, w, 0.0) - 0.5) < 1e-12  # half
    assert abs(temporal_ramp(t0 + w, t0, tfinal, w, 0.0) - 1.0) < 1e-12  # full
    assert abs(temporal_ramp(t0 + 5 * w, t0, tfinal, w, 0.0) - 1.0) < 1e-12  # interior


def test_temporal_ramp_cutoff():
    t0, tfinal, w = 0.0, 1000.0, 100.0
    assert abs(temporal_ramp(tfinal, t0, tfinal, 0.0, w) - 0.0) < 1e-12  # off at tfinal
    assert abs(temporal_ramp(tfinal - w / 2, t0, tfinal, 0.0, w) - 0.5) < 1e-12
    assert abs(temporal_ramp(tfinal - w, t0, tfinal, 0.0, w) - 1.0) < 1e-12
    assert abs(temporal_ramp(tfinal - 5 * w, t0, tfinal, 0.0, w) - 1.0) < 1e-12


def test_temporal_ramp_onset_and_cutoff_compose():
    t0, tfinal, w = 0.0, 1000.0, 100.0
    # Deep interior: both onset and cutoff are 1.
    assert abs(temporal_ramp(500.0, t0, tfinal, w, w) - 1.0) < 1e-12
    # Monotonic non-negative, bounded by 1, everywhere.
    for t in np.linspace(t0, tfinal, 50):
        r = temporal_ramp(t, t0, tfinal, w, w)
        assert 0.0 <= r <= 1.0


if __name__ == "__main__":
    test_scale_one_is_noop()
    test_scale_two_doubles_duration()
    test_anchor_invariant()
    test_scale_half_compresses_duration()
    test_model_scale_one_is_noop()
    test_model_scale_anchors_at_landfall()
    test_model_scale_direction_matches_data_storm()
    test_temporal_ramp_disabled_is_unity()
    test_temporal_ramp_onset()
    test_temporal_ramp_cutoff()
    test_temporal_ramp_onset_and_cutoff_compose()
    print("All storm_time_scale and ramp arithmetic tests passed.")
