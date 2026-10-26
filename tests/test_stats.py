import numpy as np
import pytest
from scipy import stats

from expkit import (bootstrap_diff, holm, mde, peeking_false_positive_rate, required_n_per_arm, srm_test,
                    two_prop_test)


def test_srm_balanced_and_skewed():
    assert srm_test(5000, 5000)[1] == pytest.approx(1.0)
    assert srm_test(5000, 5400)[1] < 0.001


def test_sample_size_matches_closed_form():
    # p=0.10 -> 0.12, alpha .05, power .8: (1.96+0.8416)^2*(0.09+0.1056)/0.02^2 = 3839.?
    n = required_n_per_arm(0.10, 0.02)
    assert 3830 <= n <= 3850


def test_mde_inverts_sample_size():
    n = required_n_per_arm(0.2, 0.01)
    assert mde(0.2, n) == pytest.approx(0.01, rel=0.05)


def test_two_prop_matches_statsmodels_style_z():
    r = two_prop_test(100, 1000, 130, 1000)
    pool = 230 / 2000
    z = 0.03 / np.sqrt(pool * (1 - pool) * 2 / 1000)
    assert r.z == pytest.approx(z)
    assert r.p_value == pytest.approx(2 * stats.norm.sf(z))
    assert r.ci_low < r.diff < r.ci_high


def test_newcombe_known_value():
    # Newcombe (1998) Table II example: 56/70 vs 48/80 -> diff 0.2, CI (0.0524, 0.3339)
    r = two_prop_test(48, 80, 56, 70)
    assert r.diff == pytest.approx(0.2)
    assert r.ci_low == pytest.approx(0.0524, abs=5e-4)
    assert r.ci_high == pytest.approx(0.3339, abs=5e-4)


def test_holm():
    assert holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])


def test_bootstrap_contains_true_diff():
    rng = np.random.default_rng(1)
    a, b = rng.normal(0, 1, 2000), rng.normal(0.3, 1, 2000)
    pt, lo, hi = bootstrap_diff(a, b, n_boot=500)
    assert lo < 0.3 < hi


def test_peeking_inflates_false_positives():
    fixed, peek = peeking_false_positive_rate(0.2, 2000, looks=10, sims=1500)
    assert 0.03 < fixed < 0.07
    assert peek > 2 * fixed
