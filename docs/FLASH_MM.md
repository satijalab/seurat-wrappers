Running FLASH-MM on Seurat objects for differential expression analysis
================
Compiled: May 01, 2026

- [Dataset: interferon-stimulated and control
  PBMC](#dataset-interferon-stimulated-and-control-pbmc)
- [Preprocessing: filtration and
  normalization.](#preprocessing-filtration-and-normalization)
- [Differential expression analysis using
  FLASH-MM](#differential-expression-analysis-using-flash-mm)
  - [Filtering cell types and genes](#filtering-cell-types-and-genes)
  - [Model design: fixed and random effect
    formulas](#model-design-fixed-and-random-effect-formulas)
  - [Fitting LMMs by RunFLASHMM](#fitting-lmms-by-runflashmm)
  - [Fitting LMMs with multiple random-effects
    components](#fitting-lmms-with-multiple-random-effects-components)

This vignette demonstrates how to run FLASH-MM on Seurat objects for
single-cell differential expression (DE) analysis. For demonstration
purposes, we will use the interferon-beta stimulated human PBMC dataset
(ifnb) to perform DE analysis within the same cell type across
conditions, as described in the Seurat [differential expression
testing](https://satijalab.org/seurat/articles/de_vignette) vignette.

FLASH-MM is a fast and scalable DE analysis method based on linear
mixed-effects models, which help address intra-subject correlation and
inter-subject variability in single-cell RNA-seq data. See [FLASHMM
vignette](https://cran.r-project.org/package=FLASHMM) for details. If
you find FLASH-MM useful for your publication, please cite:

> *FLASH-MM: fast and scalable single-cell differential expression
> analysis using linear mixed-effects models*
>
> Changjiang Xu, Delaram Pouyabahar, Veronique Voisin, Hamed Heydari &
> Gary D. Bader
>
> Nature Communications 17, 2026
>
> doi:
> [10.1038/s41467-026-69063-2](https://doi.org/10.1038/s41467-026-69063-2)
>
> GitHub: <https://github.com/BaderLab/FLASHMM>

Prerequisites to install:

- [Seurat](https://satijalab.org/seurat/install)
- [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
- [SeuratData](https://github.com/satijalab/seurat-data)
- [FLASHMM](https://cran.r-project.org/package=FLASHMM)

``` r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(FLASHMM)
```

## Dataset: interferon-stimulated and control PBMC

The interferon-beta stimulated human PBMCs dataset (ifnb) is available
via the [SeuratData](https://github.com/satijalab/seurat-data) package.

``` r
InstallData("ifnb")
ifnb <- LoadData("ifnb")
ifnb
#> An object of class Seurat 
#> 14053 features across 13999 samples within 1 assay 
#> Active assay: RNA (14053 features, 0 variable features)
#>  2 layers present: counts, data

# download donor (sample) information from the Github repo of the source data
info_repo <- "https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/"
ctrl <- read.table(url(paste0(info_repo, "ye1.ctrl.8.10.sm.best")), head = T, stringsAsFactors = F)
stim <- read.table(url(paste0(info_repo, "ye2.stim.8.10.sm.best")), head = T, stringsAsFactors = F)
info <- rbind(ctrl, stim)

# rename the cell IDs by substituting the '-' into '.'
info$BARCODE <- gsub("-", ".", info$BARCODE)
# only keep the cells with high-confidence sample ID (singlets, SNG)
info <- info[grep("SNG", info$BEST), ]
# remove cells with duplicated IDs in both ctrl and stim groups
info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]

# add the sample IDs to Seurat Object ifnb
rownames(info) <- info$BARCODE
ifnb <- AddMetaData(ifnb, metadata = info[, "BEST", drop = F], col.name = "donor_id")
# remove cells without donor IDs (cells that couldn't be mapped to a donor)
ifnb <- subset(ifnb, subset = !is.na(donor_id))
dim(ifnb)
#> [1] 14053 13668
```

## Preprocessing: filtration and normalization.

The standard pre-processing workflow for scRNA-seq data in Seurat,
including filtration and normalization, is described in [Seurat - Guided
Clustering
Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial).

``` r
# We filter cells that have unique feature counts over 2,500 or less than 200
ifnb <- subset(ifnb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# We normalize the data using the log-normalization method, i.e., log(1 +
# 10000*counts/library.size).
ifnb <- NormalizeData(ifnb, normalization.method = "LogNormalize", scale.factor = 10000)
```

## Differential expression analysis using FLASH-MM

In this example, we are interested in identifying what genes change in
different conditions for cells of the same type.

### Filtering cell types and genes

We filter out the cell types with a low number of cells and the genes
expressed very lowly.

``` r
# add log(library-size) to Seurat Object
ifnb <- AddMetaData(ifnb, metadata = log(colSums(ifnb[["RNA"]]$counts)), col.name = "loglib")

# rename 'seurat_annotations' as 'celltype'
colnames(ifnb@meta.data)[grep("annotations", colnames(ifnb@meta.data))] <- "celltype"

# filtering out cell types with a low number of cells
min_cells <- 100
celltype_counts <- table(ifnb$celltype)
keep_types <- names(celltype_counts)[celltype_counts >= min_cells]
ifnb <- subset(ifnb, subset = celltype %in% keep_types)

# keep genes expressed in at least 5 cells in each cell-type and treatment
celltype_trt <- paste0(ifnb@meta.data$celltype, ifnb@meta.data$stim)
gene_counts_celltypes <- do.call(cbind, tapply(1:ncol(ifnb), celltype_trt, function(i) rowSums(ifnb[["RNA"]]$counts[,
    i] > 0)))

nCells <- 5
keep_genes <- rownames(gene_counts_celltypes)[apply(gene_counts_celltypes >= nCells, 1, all)]
ifnb <- subset(ifnb, features = keep_genes)
dim(ifnb)
#> [1]  1929 13612
```

### Model design: fixed and random effect formulas

We define a fixed-effects model to capture systematic effects of
logarithm of the library size (loglib) , cell type (celltype), and the
interaction between cell type and IFN stimulation (stim), while a
random-effects model is included to account for variability between
donors (donor_id). The fixed and random effect variables in the model
formulas cannot contain special characters such as hyphen `'-'` and
empty `' '`. We rename the donor (sample) IDs by substituting `'-'` with
`''`, and the cell-types by substituting `' '` with `'.'`.

``` r
ifnb@meta.data$donor_id <- gsub("\\-", "", ifnb@meta.data$donor_id)
ifnb@meta.data$celltype <- gsub(" ", ".", ifnb@meta.data$celltype)

# fixed-effects formula
fixed <- as.formula(~0 + loglib + celltype + celltype:stim)

# random-effects formulas
random <- as.formula(~0 + donor_id)
```

### Fitting LMMs by RunFLASHMM

We run FLASHMM on the Seurat object ifnb using the fixed and random
effects we specified, with is.counts = FALSE since we are working with
normalized data.

``` r

# use normalized data
t0 <- Sys.time()
fit <- RunFLASHMM(ifnb, fixed, random, is.counts = FALSE, verbose = FALSE)
difftime(Sys.time(), t0)
#> Time difference of 12.51673 secs

# p-values
fit$p[, 1:3]
#>                                      NOC2L         ISG15         SDF4
#> loglib                        4.663435e-05  2.794627e-03 3.542139e-03
#> celltypeB                     9.721588e-03  6.181113e-01 1.094274e-01
#> celltypeB.Activated           6.017404e-01  2.338908e-01 4.309846e-01
#> celltypeCD14.Mono             3.099491e-04  8.649211e-01 1.631960e-01
#> celltypeCD16.Mono             2.909658e-03  6.539084e-02 7.596284e-02
#> celltypeCD4.Memory.T          9.603533e-03  4.799500e-01 2.797249e-01
#> celltypeCD4.Naive.T           2.456729e-02  2.806173e-01 6.960608e-02
#> celltypeCD8.T                 1.487690e-02  7.262973e-01 1.313003e-01
#> celltypeDC                    5.915536e-03  9.387579e-01 6.041559e-01
#> celltypeMk                    1.619308e-02  8.512751e-01 1.413589e-01
#> celltypeNK                    3.047048e-03  5.111534e-01 3.890581e-01
#> celltypepDC                   1.576615e-02  3.366412e-01 8.329438e-02
#> celltypeT.activated           8.207790e-02  4.278238e-01 4.331920e-01
#> celltypeB:stimSTIM            9.796296e-01  0.000000e+00 3.067832e-01
#> celltypeB.Activated:stimSTIM  5.055863e-04  0.000000e+00 6.397483e-02
#> celltypeCD14.Mono:stimSTIM    5.891962e-03  0.000000e+00 3.429434e-10
#> celltypeCD16.Mono:stimSTIM    4.083960e-03  0.000000e+00 2.558620e-01
#> celltypeCD4.Memory.T:stimSTIM 4.502681e-02  0.000000e+00 1.873564e-02
#> celltypeCD4.Naive.T:stimSTIM  2.749151e-01  0.000000e+00 1.290844e-01
#> celltypeCD8.T:stimSTIM        1.226169e-01  0.000000e+00 1.453684e-01
#> celltypeDC:stimSTIM           9.065234e-03  0.000000e+00 1.077434e-01
#> celltypeMk:stimSTIM           4.989064e-01 7.976078e-180 3.716982e-01
#> celltypeNK:stimSTIM           8.607877e-01  0.000000e+00 7.368070e-01
#> celltypepDC:stimSTIM          7.126400e-01 1.284493e-142 4.336023e-03
#> celltypeT.activated:stimSTIM  4.812262e-01  0.000000e+00 8.431043e-02
```

### Fitting LMMs with multiple random-effects components

To perform differential expression analysis between treatments (control
versus interferon-beta stimulation) across all cell types, we model
donor ID and cell type as random effects, as defined below.

``` r
# fixed-effects formula
fixed <- as.formula(~loglib + stim)

# Two random-effects components
random <- list(~0 + donor_id, ~0 + celltype)

t0 <- Sys.time()
fit <- RunFLASHMM(ifnb, fixed, random, is.counts = FALSE, verbose = FALSE)
difftime(Sys.time(), t0)
#> Time difference of 15.50337 secs

# p-values
fit$p[, 1:4]
#>                    NOC2L        ISG15         SDF4       UBE2J2
#> (Intercept) 3.755855e-02 1.988020e-05 3.298881e-01 1.615337e-02
#> loglib      1.413369e-04 7.811854e-12 3.520181e-03 9.244078e-06
#> stimSTIM    6.319996e-05 0.000000e+00 8.159833e-11 1.583305e-15
```
