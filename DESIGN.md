# Design notes

## Methods
- **SRM check:** χ² goodness of fit of the arm counts against the planned split, checked before any effect is read (Fabijan et al., 2019, *Diagnosing Sample Ratio Mismatch in Online Controlled Experiments*, KDD).
- **Sample size:**

  n = (z₁₋α/₂ + z₁₋β)² · [p₁(1−p₁) + p₂(1−p₂)] / Δ²

  The MDE uses the variance at the baseline. The test checks this against the closed form.
- **Test:** pooled two-proportion z-test.
- **Confidence interval:** Newcombe hybrid score interval (method 10; Newcombe, 1998, *Statistics in Medicine* 17:873-890). The Wald interval under-covers near 0 and 1. A unit test reproduces Newcombe's published example: 0.2, CI [0.0524, 0.3339].
- **Guardrail:** `sum_gamerounds` is extremely skewed (max 49,854; p99 493). It is tested with Mann-Whitney, plus a percentile-bootstrap CI on the mean winsorised at p99. The raw mean is dominated by one player.
- **Multiplicity:** Holm-Bonferroni (Holm, 1979) across segment tests.
- **Peeking:** vectorised A/A simulation with pooled z-tests at k equally spaced looks. It quantifies the inflation that sequential designs (Johari et al., 2017, always-valid p-values) are built to remove.

## Decisions and alternatives rejected
- **Segmenting on `retention_1`:** a post-treatment variable. It is kept only for diagnosis, explicitly labelled non-causal. Segmenting on `sum_gamerounds` was rejected: it is also post-treatment, and it is mechanically affected by the gate.
- **statsmodels:** rejected as a dependency. The formulas are short and are tested against published values, which makes the method auditable.
- **Bayesian analysis:** not used. The decision rule (ship or keep at α = 0.05) is frequentist, and adding a second framework would not change the decision.

## Changes made
- The first peeking implementation looped in Python and conflated the fixed-horizon and peeking error rates. It was rewritten as a vectorised cumulative-sum simulation.
- An early draft stated "SRM passes". The test output (p = 0.0086) contradicted it, and the decision now leads with the SRM failure.

## Limitations
- No pre-period data, so no CUPED.
- The dataset is a public extract; the assignment logs are not available.
- The power figure is a post-hoc MDE, not a pre-registered plan.
