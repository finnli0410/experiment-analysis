"""Reusable A/B test analysis: SRM, power, proportion tests with CIs, bootstrap, segments, peeking simulation."""
from .stats import (srm_test, required_n_per_arm, mde, two_prop_test, bootstrap_diff, mann_whitney,
                    segment_table, holm, peeking_false_positive_rate, PropResult)

__all__ = ["srm_test", "required_n_per_arm", "mde", "two_prop_test", "bootstrap_diff", "mann_whitney",
           "segment_table", "holm", "peeking_false_positive_rate", "PropResult"]
