# Donor-aware epithelial marker review

An exploratory analysis inspired by the **Franklin W. Huang Lab at UCSF**. It asks which epithelial expression contrasts deserve checking in additional patients before being interpreted as shared cell-state features.

This is an analysis of a published processed-count object. It retains the authors' annotations and makes no new clinical, causal, or drug-response claims. See [QA.md](QA.md) for the completed numerical checks and their limits.

## Research purpose

Song and colleagues describe heterogeneous epithelial states in localized prostate cancer and incomplete representation of those states in organoid cultures. A useful next question is whether a selected marker contrast is shared across donors or disproportionately influenced by donors contributing many cells. This is an inferred research opportunity, not a statement about the lab's internal priorities.

The workflow compares Club cells with luminal epithelial (LE) or basal epithelial (BE) cells, or LE with BE, within a selected specimen stratum. It shows donor coverage, cell-count-weighted and equal-donor contrasts, and sensitivity to omitting one donor. Its output is a transparent marker-review list, not a diagnostic classifier or a drug-target recommendation.

Tissue–organoid retention is **not** evaluated unless verified parent-donor metadata and a corresponding analysis are added. Finding similar marker expression alone would not establish model fidelity, lineage relationships, or causal state transitions.

## Sources and provenance

- Song H, Weinstein HNW, Allegakoen P, and colleagues. *Single-cell analysis of human primary prostate cancer reveals the heterogeneity of tumor-associated epithelial cell states*. Nature Communications **13**, 141 (2022). [Article and supplements](https://doi.org/10.1038/s41467-021-27322-4), [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8748675/).
- [GEO GSE176031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176031): study sequencing data and expression matrices. This demo is intended to analyze a lab-linked processed object, not to repeat FASTQ processing.
- [Lab repository](https://github.com/franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples), inspected at commit [`71dff9b7dc3d096b6ea2da1ed2cc5bae22876982`](https://github.com/franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples/tree/71dff9b7dc3d096b6ea2da1ed2cc5bae22876982). The paper links the corresponding [author repository](https://github.com/angelussong/scRNA-seq-Analysis-of-Prostate-Cancer-Samples).
- The pinned repository [README](https://github.com/franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples/blob/71dff9b7dc3d096b6ea2da1ed2cc5bae22876982/README.md) links the [processed Seurat-object archive](https://www.dropbox.com/sh/ga770h1ik75l1o2/AADfmpEsHp8OLTtfkZaCxPd8a?dl=0), describes `active.ident` as the annotation and `orig.ident` as the patient identifier, and notes shortened cell names. These descriptions must be checked against the actual object; identifiers are not silently repaired.
- Candidate markers are author-chosen inputs from [`Seqwell_Combined_E.R`](https://github.com/franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples/blob/71dff9b7dc3d096b6ea2da1ed2cc5bae22876982/Seqwell_Combined_E.R): the epithelial marker lists around lines 232/252 and explicit club-comparison labels around lines 836–847. The fixed union was selected before inspecting effect estimates. These are not novel marker discoveries or an unbiased transcriptome-wide target screen.
- The related [Clubcell_Organoid repository](https://github.com/franklinhuanglab/Clubcell_Organoid/tree/ac2f6582ce502b33245aa35e91c7f9c1ac9b8590) was inspected at commit `ac2f6582ce502b33245aa35e91c7f9c1ac9b8590`. Inspection does not imply that its full analysis was run or that organoid data are included here.

| Release provenance item | Status |
| --- | --- |
| Selected source | `dge_E.rds`, from the repository-linked [processed-object archive](https://www.dropbox.com/sh/ga770h1ik75l1o2/AADfmpEsHp8OLTtfkZaCxPd8a?dl=1) |
| RDS SHA-256 | `c0cbde35277b2f250ee13160afe9d163b62c720dca8c9f43f121a47e98b29430` |
| Downloaded archive SHA-256 | `c969b8903d2387d5f186b9bf233cab41da2ff7f29eaea970365efed75ca5a790` |
| Assay and source metadata | Deposited `RNA` counts/data; `orig.ident` source identifiers; `sampletype` tissue stratum; `ID` state labels agree with `active.ident` |
| Source coverage | 13,322 epithelial cells, 21,877 RNA features, 15 source groups mapped to 11 donors, 73 donor/tissue/state groups |
| Panel coverage | All 20 prespecified genes present; original `IL8` is present, so no `CXCL8` alias substitution was made |
| Normalization check | Recomputed selected-gene values agree with the stored normalized assay to maximum absolute difference `5.55e-17`; whole-assay count totals agree exactly with `nCount_RNA` |
| Reproduction implementation | [Processed-object extractor](analysis/extract_patient_markers.py), [descriptive comparison script](analysis/summarize_markers.py), and [pinned dependencies](analysis/requirements.txt); no stored R commands are evaluated |
| Derived JSON SHA-256 | `783a9788e5cc9f1afa9fac96fe8397c5266ea1266a1b056af72ef52b52d04abe` |
| Verification | [Source manifest](analysis/source_manifest.json), [numerical QA](QA.md), [10 Python regression tests](analysis/test_analysis.py), and [browser-calculation tests](scripts/test_analysis.cjs); implementation checks are not biological validation |

The fixed marker panel is `KRT5`, `KRT15`, `KRT17`, `TP63`, `KLK3`, `KRT8`, `KRT18`, `KRT13`, `APOBEC3A`, `CSTB`, `LYPD3`, `SERPINB1`, `KLK2`, `PIGR`, `MMP7`, `CP`, `SCGB1A1`, `CD74`, `LTF`, and `IL8`.

Source labels comprise 5,647 LE, 2,408 ERG-negative tumor, 2,250 ERG-positive tumor, 2,145 BE, and 872 Club cells. These are the authors' labels, not new malignancy assignments. The tumor-labeled groups remain in source coverage metadata but are not the selectable state contrasts in this prototype.

### Donor and specimen mapping

`PR5249_N/T`, `PR5251_N/T`, `PR5254_N/T`, and `PR5261_N/T` map to their four respective donor roots. `PR5186`, `PR5196`, `PR5199`, and `PR5269` remain four separate radical-prostatectomy donors. `AUG_PB1`, `MAY_PB1`, and `MAY_PB2` remain three separate biopsy donors. Original identifiers are preserved; the later upstream script recode of `PR5261` to `PR5254` is not applied.

The 15 `orig.ident` values are **source sample groups, not 15 physical specimens or 15 independent patients**. The publication describes six biopsy specimens from three patients, while the deposited object presents those biopsy cells under three patient-level source groups. This object cannot recover biopsy-level variation that has already been collapsed; biopsy specimens are not counted as separate biological replicates.

The interface offers three specimen strata:

- **Radical-prostatectomy tumor:** source `Paired_T` (3,691 cells) and `Unpaired_T` (3,193 cells), combined as one sampling context while retaining source provenance.
- **Biopsy:** source `Biopsy` (3,018 cells), analyzed separately from prostatectomy specimens.
- **Grossly normal prostatectomy tissue:** source `Paired_N` (3,420 cells), not labeled healthy-control tissue.

These are source-wide counts before the chosen state pair and coverage filter. The interface and export report the smaller eligible comparison counts and all omitted donors. A donor contributes at most one paired state contrast to any selected stratum.

## Analysis contract

1. Keep source identifiers, tissue strata, donor assignments, and cell-state annotations. Multiple samples or biopsies from one donor do not become independent donors.
2. Require both comparison states within each eligible donor and selected specimen stratum. The default minimum is 20 cells per state, with 10/20/30/50-cell sensitivity settings. These are exploratory coverage choices, not established statistical sufficiency thresholds.
3. Recompute per-cell `ln(1 + 10000 × gene count / total counts across all 21877 deposited RNA features)`, then average within each donor/state group. A contrast of these means is a **mean-log-expression difference**, not a log fold change. The denominator is the full deposited RNA assay, not the 20-gene panel.
4. Compute detection as the fraction of included cells with count greater than zero. Missing genes or unavailable groups remain unavailable, not zero.
5. Compare pooled cell-count-weighted estimates with the mean of paired donor effects using exactly the same eligible donors and cells. Show when unequal cell contribution changes the interpretation.
6. Omit each eligible donor in turn and recalculate the equal-donor estimate. This tests sensitivity to donor omission; it is not held-out prediction, a confidence interval, or independent replication. Do not promote a marker into the stability shortlist with fewer than four eligible donors.
7. Export the selected settings, coverage, estimates, and exclusions alongside the review list. Do not hide missing comparison groups or rank missing values as weak biological effects.

The positive-contrast candidate rule requires at least four eligible donors, a strictly positive paired difference in at least 75% of them, a strictly positive mean after leaving out each donor, and an equal-donor mean difference meeting the selected minimum (0.10, 0.25, or 0.50). These are descriptive review settings, not validated biomarker criteria. Negative contrasts remain visible rather than being discarded.

For the default tumor-RP, Club-minus-LE comparison at 20 cells per state, four donors contribute 268 Club and 1,607 LE cells. Eight of the 20 prespecified markers meet the exploratory 0.25-threshold rule. The unrounded value controls eligibility: KRT5's mean difference is approximately 0.249667, so it does not pass even though a two-decimal display reads 0.25. [QA.md](QA.md) records baseline results and edge-case checks.

## Reproduce the analysis

The published summaries can be checked without downloading the RDS object or installing extraction dependencies:

```bash
python3 analysis/test_analysis.py
node scripts/test_analysis.cjs
python3 analysis/summarize_markers.py --data patient-marker-data.json --output /tmp/huang-marker-summary.json
```

The Python suite contains ten tests and covers all 108 supported setting combinations: three specimen policies, three state contrasts, four cell-count settings, and three effect thresholds. The browser test checks independently calculated baseline results, 12 edge cases, CSV contents/quoting, and equality between the JSON and browser bundle. The full independent comparison checked 2,160 gene–setting combinations; details and the distinction between reusable tests and source-RDS verification are in [QA.md](QA.md).

To rebuild summaries from the deposited object, use Python 3.12 and an isolated environment. The pinned extraction environment is `rdata==1.1.0`, `numpy==2.5.4`, `scipy==1.18.1`, `pandas==3.0.6`, and `xarray==2026.9.0`. The descriptive summary and bundle scripts use the Python standard library.

Obtain `dge_E.rds` from the original repository-linked archive. Keep it outside the public site; it is not bundled here. The extractor verifies its expected SHA-256, converts the relevant metadata and sparse RNA assay, and recomputes selected-gene summaries without executing stored R commands. It does not repeat alignment, original filtering, clustering, or cell annotation.

From this repository's root:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r analysis/requirements.txt
.venv/bin/python analysis/extract_patient_markers.py --rds /path/to/dge_E.rds --output patient-marker-data.json
.venv/bin/python analysis/bundle_data.py --data patient-marker-data.json --output data.js
.venv/bin/python analysis/summarize_markers.py --data patient-marker-data.json --output /tmp/huang-marker-summary.json
.venv/bin/python analysis/test_analysis.py
node scripts/test_analysis.cjs
```

`/path/to/dge_E.rds` is an explicit placeholder for the downloaded source object, not a bundled file. `patient-marker-data.json` stores the derived donor/state summaries and provenance. The comparison script writes descriptive results to the requested output path. The bundle step wraps the JSON as `window.HUANG_DATA` in `data.js` without changing its values; changing the JSON alone does not refresh that separate bundle.

To view the static page locally, run `python3 -m http.server 8000` from the repository root and open `http://localhost:8000/`. Keep downloaded source objects, virtual environments, and cell-level intermediate files outside any public deployment.

## Reuse with another cohort

For reuse, prepare the following long-table schema, or map it explicitly to the bundled compact JSON format. A public upload/parser for arbitrary new cohorts is not claimed.

| Field | Meaning |
| --- | --- |
| `gene` | Recorded gene identifier; retain an explicit alias map |
| `donor` | Independent biological donor identifier |
| `source_sample` | Original source-group identifier; record true specimen/library replicates separately when available, rather than treating this study's collapsed groups as individual specimens |
| `tissue` | Verified tissue stratum or collection context |
| `state` | Recorded cell-state annotation with its annotation method |
| `n_cells` | Integer number of cells contributing to this group |
| `mean_log_expression` | Mean expression on the explicitly documented scale |
| `detection_fraction` | Fraction of those cells with a nonzero source count |

Supply normalization and filtering metadata separately: assay name, library-size denominator, transform, scale factor, gene annotation version, state-definition method, and any batch correction. Reject duplicate keys, nonfinite values, invalid counts/fractions, inconsistent donor mappings, and missing comparison groups. If a donor has multiple samples in one stratum, document aggregation before comparing donors and do not average sample means without their cell counts. Do not merge studies with incompatible normalization under a single expression scale.

The reusable outputs are the settings-specific marker-review CSV, per-donor evidence, and a list of omitted groups with reasons. The public interface does not imply an implemented user-upload feature. The extractor and regression checks are intentionally pinned to this published object: adapting them to a new cohort requires an explicit mapping/normalization review, new provenance, and new cohort-specific expected values. Do not bypass checksum or coverage tests to make an incompatible input appear validated.

## A validation experiment that distinguishes the explanations

**Shared-state explanation:** the selected marker difference has the same direction in additional independent donors after using a state definition not based on that marker. **Donor-bias explanation:** the pooled contrast depends on a small subset of donors, collection conditions, or cell composition and does not replicate consistently.

Preselect both apparently stable and donor-sensitive markers before collecting the validation cohort. In independent specimens, obtain both comparison states from each donor using an orthogonal annotation and measure the selected markers with an appropriate validated assay, such as targeted RNA measurement or protein staining when a validated reagent is available. Include positive and negative markers, assay controls, recorded collection/processing conditions, and balanced processing batches. Evaluate the paired donor effects, not thousands of cells as independent replicates.

Consistent within-donor replication would support further study of the marker as a state-associated feature. Disappearance or reversal across independent donors would downgrade it. Disagreement between RNA and protein would motivate assay and biological follow-up rather than a claim of success. This experiment remains a proposal; it has not been performed by this demo.

## Interpretation limits

- Published state annotations partly use marker expression, so contrasts of those markers are partly circular. This analysis reviews consistency; it does not independently validate the states.
- Cells are not independent biological replicates. Coverage filters, equal-donor weighting, and leave-one-donor-out checks do not remove all sampling, batch, or annotation bias.
- The source paper's tumor identification partly relies on inferred copy-number changes. Shallow sequencing can miss focal alterations, and grossly normal specimens may contain tumor cells.
- If external healthy-reference data are used in a future extension, study and processing differences remain potential confounders. This draft does not claim such a comparison was performed.
- Expression or detection does not establish target essentiality, drug response, prognosis, lineage, or causality. No cell-level significance tests, clinical predictions, or therapeutic recommendations are proposed here.
- This is a focused analysis of processed data, not reproduction of the full upstream R workflow or the later prostate atlas.

## Contribution and permissions

Cathy's contribution is framed around data curation, interpretation of single-cell contrasts, and transparent follow-up reasoning. The computational implementation is assistant-supported. It does not establish unaided programming proficiency, independent organoid culture, or new hands-on experimental experience.

Demo code is distributed under [GNU GPL version 3](LICENSE). Upstream attribution and separate article/data terms are described in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). No full source RDS objects, count matrices, or raw sequencing files are redistributed in this repository.
