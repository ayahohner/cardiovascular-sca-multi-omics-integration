import numpy as np

from cardio_sca_mvp.signals import benjamini_hochberg


def test_benjamini_hochberg_is_monotone_in_p_value_order() -> None:
    p_values = np.array([0.01, 0.04, 0.03, 0.002, np.nan])
    adjusted = benjamini_hochberg(p_values)
    finite_order = np.argsort(p_values[:-1])
    assert np.all(np.diff(adjusted[:-1][finite_order]) >= 0)
    assert np.isnan(adjusted[-1])
    assert np.all((adjusted[:-1] >= p_values[:-1]) & (adjusted[:-1] <= 1.0))
