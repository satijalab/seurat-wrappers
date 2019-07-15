# SeuratWrappers

SeuratWrappers is a collection of community-provided methods and extensions for [Seurat](https://satijalab.org/seurat/), curated by the Satija Lab at NYGC. These methods comprise functionality not presently found in Seurat, and are able to be updated much more frequently.

Please see our [contribution guide](https://github.com/satijalab/seurat.wrappers/wiki) for assistance and guidelines in developing and adding new methods to SeuratWrappers

Individual method vignettes can be found in the [`docs/`](https://github.com/satijalab/seurat.wrappers/tree/master/docs) directory, we recommend looking at the standard markdown (`*.md`) files when viewing on GitHub

Installation can be accomplished through [devtools](https://cran.r-project.org/package=devtools)

```R
devtools::install_github('satijalab/seurat.wrappers')
```

## Method Listing

| Package | Vignette | Reference | Source |
| ------- | -------- | --------- | ------ |
| Conos | [Integration of datasets using Conos](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/conos.html) | Barkas et al, Nature Methods 2019 | https://github.com/hms-dbmi/conos |
| LIGER | [Integrating Seurat objects using LIGER](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/liger.html) | Welch et al, Cell 2019 | https://github.com/MacoskoLab/liger |
| fastMNN | [Running fastMNN on Seurat Objects](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.html) | Nature Biotechnology 2018 | https://bioconductor.org/packages/release/bioc/html/scran.html |
| Harmony | [Integration of datasets using Harmony](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/harmony.html) | Korsunsky et al, bioRxiv 2018 | https://github.com/immunogenomics/harmony |
| ALRA | [Zero-preserving imputation with ALRA](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/alra.html) | Linderman et al, bioRxiv 2018 | https://github.com/KlugerLab/ALRA |
| Velocity | [Estimating RNA Velocity using Seurat](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/velocity.html) | La Manno et al, Nature 2018 | https://velocyto.org |
