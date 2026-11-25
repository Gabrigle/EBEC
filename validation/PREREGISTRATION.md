# EBEC-scan: pre-specified validation criteria

Fixed on 2026-10-19, BEFORE running any simulation of the new statistic.
Do not edit after seeing results. If a criterion must change, add a dated
amendment at the end and report both versions.

## Statistic under evaluation
`ebec_scan(y, p, q = c(0.01, 0.05, 0.10, 0.20), classes = TRUE)`:
max standardized excess over total loss, tail sums of the k = ceil(q n) worst
losses, and the same tail sums for false negatives (events) and false positives
(non-events). Comparators: Spiegelhalter, Hosmer-Lemeshow (10 groups), Cox
calibration LR test (intercept and slope, 2 df), and the same scan without
class components (`scan_tot`). alpha = 0.05.

## Criteria
- **C1 (size, blocking).** On H0 cells, the pooled rejection rate of `scan_all`
  is at most 0.05 + 3 SE (SE = sqrt(0.05 * 0.95 / N)) and at least 0.03, and at
  most 10% of cells have a Wilson 95% interval entirely above 0.05.
- **C2 (rare confident errors, blocking for the "new test" claim).** In the 12
  rare-confident-error scenarios, `scan_all` power exceeds Spiegelhalter power by
  at least 0.15 in at least 8 of 12 scenarios (power averaged over the cells run).
- **C3 (class components, blocking for any FN/FP claim).** In the 8 scenarios
  that target false negatives (4 FN under-prediction + 4 rare FN-only) and in the
  8 that target false positives, `scan_all` power exceeds `scan_tot` power by at
  least 0.05 in at least half of each group.
- **C4 (global shifts, reported, not blocking).** Power of `scan_all` minus
  Spiegelhalter in the 2 global-overconfidence scenarios, reported as is.
- **C5 (fitted models, reported).** Type I error of `ebec_scan` on in-sample fitted
  probabilities versus 5-fold cross-validated probabilities (amended, see below).

## Decision rule
- C1 fails: stop; the test is not valid; fix before anything else.
- C1 passes, C2 and C3 pass: claim a valid test with demonstrated advantage in the
  rare-confident-error regime, and class-specific components.
- C1 passes, C2 passes, C3 fails: claim the advantage in the rare-error regime;
  do NOT claim FN/FP attribution.
- C1 passes, C2 fails: no claim of added detection power; present as a valid
  descriptive diagnostic of tail loss.

## Scenarios, all reported (no selection)
H0: n in {100, 250, 500, 1000} x prevalence in {0.10, 0.30, 0.50}.
H1: the 14 scenarios of the 18-Oct run, plus 12 rare-confident-error scenarios
(c in {0.01, 0.02}, stated probability extreme = 0.02 or 0.10 from the boundary,
type in {mixed, FN-only, FP-only}), n in {250, 500, 1000}, same prevalences.
In every scenario the stated probabilities p are generated without using y.

## Amendment 1 (2026-10-19, after the pilot run, before any full run)
The pilot showed that `ebec_scan_glm` (parametric-bootstrap refit) is degenerate:
the statistic depends almost only on the fitted probabilities, so the refit null
is far too wide (type I 0/200, median p = 0.89). The function was removed from the
package. C5 now compares in-sample probabilities with 5-fold cross-validated
probabilities. C1-C4 are unchanged.

## Amendment 2 (2026-10-20, after the full run of v2, before any run of v3)
Result of the full v2 run (248 cells, 0 errors): C1 passed (scan_all 0.046), C2
passed (10 of 12), C3 failed: the false-negative components gave no gain over
`scan_tot` (FN-targeting scenarios 0 of 8 with a gain >= 0.05; e.g. FN
under-prediction c = 0.10, d = 3: 0.915 vs 0.918), while the false-positive
components did (FP-targeting scenarios: 1.00).
Decision, taken from the v2 results and therefore to be re-validated with new seeds:
1. Default statistic: `ebec_scan(..., classes = "fp")`: total loss, tail sums, and
   FP tail sums. FN components are available only with `classes = "both"`
   (exploratory); no FN-specific claim is made.
2. The v3 validation (`EBEC3_validation.R`, seed 20261020, same grid) evaluates
   `scan_fp` as the primary statistic; `scan_both` and `scan_tot` are reported.
3. Criteria for v3: C1 and C2 as before, with `scan_fp` in place of `scan_all`.
   C3 becomes FP-only: in the 8 FP-targeting scenarios (4 FP over-prediction + 4
   rare FP-only), `scan_fp` exceeds `scan_tot` by >= 0.05 in at least half. No
   criterion for FN. C4 and C5 unchanged (reported). New C6 (reported): mean of
   `scan_fp - scan_tot` over the non-FP scenarios, tolerated >= -0.03.
4. Decision rule: C1 fails, stop. C1 and C2 pass and C3 passes: claim a valid test
   with demonstrated advantage in the rare-confident-error regime and a useful
   false-positive component; no FN claim in any case. C3 fails: claim only the
   rare-error advantage. C2 fails: valid descriptive diagnostic only.
5. v2 results stay reported as the exploratory phase that led to this amendment.

### Amendment 2b (2026-10-20, same day, before any v3 run): naming and reporting
- The statistic is exported as `ebec()` (class `ebec`). The earlier share-of-loss
  index is kept only as `ebec_share()`, `ebec_share_fn()` and `ebec_share_fp()`
  (legacy, descriptive). The computation of the statistic is unchanged; the v3
  scripts call `ebec()` with `ci = FALSE`.
- New descriptive output, not part of any criterion: ratio observed / expected and
  percentage excess per component, with a bootstrap 95% interval (cases
  resampled; expected value re-simulated from the resampled probabilities). A
  coverage check in a small fixed-design simulation (n = 400, about 8% events,
  120 replicates, three designs) gave 0.88-0.93 under H0 and a rare
  false-negative design, and 0.97-1.00 under global overconfidence. The interval
  is documented as descriptive and not adjusted for the choice of component.
