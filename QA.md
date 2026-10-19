# Numerical verification

Checked on 19 October 2026. These checks concern implementation correctness and source consistency, not biological validation, clinical performance, or endorsement by the source authors.

## Independent source checks

The original `dge_E.rds` object was read separately from the transformation script. Its SHA-256 was `c0cbde35277b2f250ee13160afe9d163b62c720dca8c9f43f121a47e98b29430`.

- All 266,440 selected marker-by-cell raw count values matched the extraction cache.
- All 13,322 library totals, summed directly from the deposited sparse-count vectors, matched the cache and metadata totals.
- All 266,440 selected per-cell normalized values were checked using `log1p(10000 × count / deposited RNA total)`.
- All 2,920 group mean-expression and detection fields agreed with a second calculation: 73 donor–specimen–state groups × 20 genes × two summaries.
- Three positive source-assay spot checks used PIGR, KLK3 and KRT5. Recomputed values matched the deposited normalized RNA assay. These checks did not rerun alignment, clustering or annotation.

## Browser calculation checks

A separately written calculation generated expectations for all 108 combinations of three specimen policies, three state contrasts, four minimum-cell settings and three effect thresholds. The browser analysis agreed across all 2,160 gene–setting combinations and 12,962 numerical/null comparisons, including weighting, detection differences, donor eligibility, leave-one-out bounds and candidate designation.

Twelve synthetic unit-test cases additionally checked empty coverage, one or three eligible donors, missing gene values, exact and near-rounded thresholds, the 75% positive-donor boundary, leave-one-out reversal, unequal cell counts, and invalid gene/policy requests. Synthetic inputs are used only for tests; they are not demo research observations.

The default tumor-RP, Club-minus-LE comparison at 20 cells per state uses four eligible patients and 1,875 cells: 268 Club and 1,607 LE cells. Eight of the fixed 20 markers meet the exploratory rules at a 0.25 mean-log-expression threshold. KRT5 does **not** pass: its unrounded difference is about 0.249667, despite rounding to 0.25 at two decimal places.

## Reusable checks

With Node.js available, run from the repository root:

```sh
node scripts/test_analysis.cjs
```

This public test checks four independently calculated baseline gene results, the candidate list, twelve edge cases, CSV values and quoting, and equality between the browser data bundle and JSON. The full source-RDS cross-check requires obtaining the documented source object; this repository does not redistribute the source object or cell-level cache. See the README for extraction and analysis commands.

Correct arithmetic cannot resolve source annotation circularity, small donor coverage, acquisition/batch differences, or the absence of an independent validation cohort. Candidate labels are descriptive follow-up priorities, not significance results or validated biomarkers.
