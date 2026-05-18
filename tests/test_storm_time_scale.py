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


if __name__ == "__main__":
    test_scale_one_is_noop()
    test_scale_two_doubles_duration()
    test_anchor_invariant()
    test_scale_half_compresses_duration()
    print("All storm_time_scale arithmetic tests passed.")
