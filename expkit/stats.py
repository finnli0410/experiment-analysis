from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats


def srm_test(n_a: int, n_b: int, expected_share_a: float = 0.5) -> tuple[float, float]:
    """Sample-ratio mismatch: chi-square goodness of fit against the planned split. Returns (chi2, p)."""
    n = n_a + n_b
    chi2, p = stats.chisquare([n_a, n_b], [n * expected_share_a, n * (1 - expected_share_a)])
    return float(chi2), float(p)


def required_n_per_arm(p_base: float, mde_abs: float, alpha: float = 0.05, power: float = 0.8) -> int:
    """Two-sided two-proportion z-test, equal arms (unpooled normal approximation)."""
    p2 = p_base + mde_abs
    za, zb = stats.norm.ppf(1 - alpha / 2), stats.norm.ppf(power)
    var = p_base * (1 - p_base) + p2 * (1 - p2)
    return int(np.ceil((za + zb) ** 2 * var / mde_abs ** 2))


def mde(p_base: float, n_per_arm: int, alpha: float = 0.05, power: float = 0.8) -> float:
    """Minimum detectable absolute difference for a given n (variance evaluated at p_base)."""
    za, zb = stats.norm.ppf(1 - alpha / 2), stats.norm.ppf(power)
    return float((za + zb) * np.sqrt(2 * p_base * (1 - p_base) / n_per_arm))


@dataclass(frozen=True)
class PropResult:
    p_a: float
    p_b: float
    diff: float          # p_b - p_a
    ci_low: float
    ci_high: float
    z: float
    p_value: float
    rel_lift: float


def two_prop_test(x_a: int, n_a: int, x_b: int, n_b: int, alpha: float = 0.05) -> PropResult:
    """Pooled z-test for H0: p_a = p_b; Newcombe (Wilson score, method 10) CI for p_b - p_a."""
    pa, pb = x_a / n_a, x_b / n_b
    pool = (x_a + x_b) / (n_a + n_b)
    se0 = np.sqrt(pool * (1 - pool) * (1 / n_a + 1 / n_b))
    z = (pb - pa) / se0
    p = 2 * stats.norm.sf(abs(z))
    zc = stats.norm.ppf(1 - alpha / 2)

    def wilson(x: int, n: int) -> tuple[float, float]:
        ph = x / n
        c = (ph + zc**2 / (2 * n)) / (1 + zc**2 / n)
        h = zc * np.sqrt(ph * (1 - ph) / n + zc**2 / (4 * n**2)) / (1 + zc**2 / n)
        return c - h, c + h

    la, ua = wilson(x_a, n_a); lb, ub = wilson(x_b, n_b)
    d = pb - pa
    lo = d - np.sqrt((pb - lb) ** 2 + (ua - pa) ** 2)
    hi = d + np.sqrt((ub - pb) ** 2 + (pa - la) ** 2)
    return PropResult(pa, pb, d, float(lo), float(hi), float(z), float(p), d / pa)


def bootstrap_diff(a: np.ndarray, b: np.ndarray, stat=np.mean, n_boot: int = 2000, seed: int = 0,
                   alpha: float = 0.05) -> tuple[float, float, float]:
    """Percentile bootstrap CI for stat(b) - stat(a). Returns (point, low, high)."""
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a), np.asarray(b)
    boots = np.array([stat(rng.choice(b, b.size)) - stat(rng.choice(a, a.size)) for _ in range(n_boot)])
    return float(stat(b) - stat(a)), float(np.quantile(boots, alpha / 2)), float(np.quantile(boots, 1 - alpha / 2))


def mann_whitney(a, b) -> float:
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def holm(pvals: list[float]) -> list[float]:
    """Holm-Bonferroni adjusted p-values (family-wise error control for segment tests)."""
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * pvals[i])
        adj[i] = min(1.0, running)
    return adj.tolist()


def segment_table(df: pd.DataFrame, arm: str, metric: str, segment: str, control, treatment) -> pd.DataFrame:
    rows = []
    for seg, g in df.groupby(segment, observed=True):
        a, b = g[g[arm] == control][metric], g[g[arm] == treatment][metric]
        r = two_prop_test(int(a.sum()), len(a), int(b.sum()), len(b))
        rows.append({segment: seg, "n_control": len(a), "n_treatment": len(b), "p_control": r.p_a,
                     "p_treatment": r.p_b, "diff_pp": r.diff * 100, "ci_low_pp": r.ci_low * 100,
                     "ci_high_pp": r.ci_high * 100, "p_value": r.p_value})
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out.p_value.tolist())
    return out


def peeking_false_positive_rate(p: float, n_per_arm: int, looks: int, sims: int = 2000, alpha: float = 0.05,
                                seed: int = 0) -> tuple[float, float]:
    """A/A simulation (both arms share true rate p). Pooled z-test at `looks` equally spaced checkpoints.
    Returns (FPR when testing only at the final look, FPR when stopping at the first significant look)."""
    rng = np.random.default_rng(seed)
    n = np.linspace(n_per_arm / looks, n_per_arm, looks).astype(int)
    zc = stats.norm.ppf(1 - alpha / 2)
    ca = np.cumsum(rng.random((sims, n_per_arm)) < p, axis=1)[:, n - 1]
    cb = np.cumsum(rng.random((sims, n_per_arm)) < p, axis=1)[:, n - 1]
    pool = (ca + cb) / (2 * n)
    se = np.sqrt(pool * (1 - pool) * 2 / n)
    z = np.divide(np.abs(cb - ca) / n, se, out=np.zeros_like(se), where=se > 0)
    sig = z > zc
    return float(sig[:, -1].mean()), float(sig.any(axis=1).mean())
