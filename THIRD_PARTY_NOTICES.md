# Third-party sources and notices

This notice distinguishes software, article content, and linked data; it does not assign a blanket license to all upstream material. Exact files, hashes, transformations, and limits are recorded in the [source manifest](analysis/source_manifest.json), [README](README.md), and [numerical QA](QA.md).

## Demo software

The demo's code is made available under GNU General Public License version 3; see [LICENSE](LICENSE). Newly written transformation and interface code must be identified as this project's implementation, not as code authored or endorsed by the Huang Lab. No full upstream R script is intended to be bundled.

## Upstream analysis software and marker inputs

**Hanbing Song and the Franklin Huang Lab / study authors**

- [scRNA-seq-Analysis-of-Prostate-Cancer-Samples](https://github.com/franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples), commit `71dff9b7dc3d096b6ea2da1ed2cc5bae22876982`.
- [Clubcell_Organoid](https://github.com/franklinhuanglab/Clubcell_Organoid), commit `ac2f6582ce502b33245aa35e91c7f9c1ac9b8590`.

Both inspected repositories contain a GNU GPL version 3 license. Retain their attribution when reusing or adapting software. The 2022 study repository supplies the processed-object link and the author-selected marker lists used to define the fixed panel. Marker names and author selections are research inputs, not discoveries made by this demo. The pinned files, source hashes, and actual usage are documented in the README. The related organoid repository was inspected for context; it is not a source of the tissue-cell observations analyzed here.

The later `scRNA_PCa_Atlas` repository was inspected during discovery but is not an input to this demo. No license was identified during that inspection; no code or dataset from it is redistributed, and no reproduction of that atlas is claimed.

## Article and supplementary materials

Song H, Weinstein HNW, Allegakoen P, and colleagues. *Single-cell analysis of human primary prostate cancer reveals the heterogeneity of tumor-associated epithelial cell states*. Nature Communications 13, 141 (2022). DOI: [10.1038/s41467-021-27322-4](https://doi.org/10.1038/s41467-021-27322-4).

The [article's license statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC8748675/) specifies [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/), subject to separately credited third-party material. Any reused article or supplementary material must credit the authors and source, link the license, and identify changes. The present narrative paraphrases the paper and adds a separately labeled exploratory analysis and validation proposal; it is not an author-approved interpretation.

## Linked processed data

The study's repository README links processed Seurat objects through Dropbox; the study also deposits sequencing data and matrices as [GSE176031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176031). This analysis uses the linked `dge_E.rds` object, SHA-256 `c0cbde35277b2f250ee13160afe9d163b62c720dca8c9f43f121a47e98b29430`.

Source identifiers are retained as deposited. The object's 15 source groups are not a count of physical specimens: six published biopsy specimens appear as three patient-level source groups. No additional specimen-level resolution or donor identities are inferred.

**No explicit, separate license for the linked RDS objects has been established.** Public accessibility does not by itself establish that the repository's software license or the article's CC BY license covers those separate files. This project does not redistribute the source RDS archive, whole expression matrices, cell-level files, or raw reads. Readers who want to reproduce the analysis should obtain the exact documented source from its original provider and review its terms.

Only minimal newly computed, non-cell-level summary statistics and provenance needed for the demo are published. These summaries remain attributed to the source study. This notice does not relabel the original data as GPL or claim unrestricted upstream redistribution rights. If later provenance checks reveal additional restrictions, the release contents must be revised accordingly.

## No endorsement

This independent exploratory project is not affiliated with or endorsed by the source authors, Franklin Huang Lab, UCSF, or the repositories cited above. Its proposed marker review and follow-up experiment are not clinical or experimentally validated recommendations.
