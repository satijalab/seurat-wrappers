Running BANKSY with Seurat
================
Compiled: April 04, 2024

- [Introduction](#introduction)
- [Overview](#overview)
- [Running BANKSY within Seurat’s spatial
  framework](#running-banksy-within-seurats-spatial-framework)
- [Running BANKSY with locations provided
  explicitly](#running-banksy-with-locations-provided-explicitly)
- [Multi-sample analysis](#multi-sample-analysis)
- [Spatial data integration with
  Harmony](#spatial-data-integration-with-harmony)
- [Getting help](#getting-help)

## Introduction

In this vignette, we describe how to run BANKSY with Seurat objects. If
you use BANKSY in your research, please cite

> *BANKSY unifies cell typing and tissue domain segmentation for
> scalable spatial omics data analysis*
>
> Vipul Singhal, Nigel Chou, Joseph Lee, Yifei Yue, Jinyue Liu, Wan Kee
> Chock, Li Lin, Yun-Ching Chang, Erica Mei Ling Teo, Jonathan Aow, Hwee
> Kuan Lee, Kok Hao Chen & Shyam Prabhakar
>
> Nature Genetics, 2024
>
> doi:
> [10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)
>
> Website: <https://prabhakarlab.github.io/Banksy>

BANKSY is a method that incorporates neighborhood information for
clustering spatial omics data. By doing so, BANKSY is able to

- improve cell-type assignment in noisy data
- distinguish subtly different cell-types stratified by microenvironment
- identify spatial domains sharing the same microenvironment

The amount of neighborhood information incorporated is controlled by a
parameter `lambda` in \[0,1\], with higher values giving more weight to
the neighbourhood information during clustering.

## Overview

The `RunBanksy` function implemented with the *SeuratWrappers* package
allows users to run BANKSY with Seurat objects. We describe two options
of running `RunBanksy`. The first is within Seurat’s spatial framework
(see [here](https://satijalab.org/seurat/articles/spatial_vignette.html)
and
[here](https://satijalab.org/seurat/articles/spatial_vignette_2.html))
and requires a Seurat object and a lambda parameter as mandatory input.
The second option works with Seurat objects that do not have spatial
information stored within, and therefore requires an additional argument
giving the locations of the cell centroids or spots.

**Caveat**: `ScaleData` should not be run after a call to `RunBanksy`;
`RunBanksy` populates the `scale.data` slot with the scaled BANKSY
matrix. Calling `ScaleData` after `RunBanksy` performs gene-wise
z-scaling, negating the effect of `lambda`.

Prerequisites to install:

- [Seurat](https://satijalab.org/seurat/install)
- [SeuratData](https://github.com/satijalab/seurat-data)
- [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
- [Banksy](https://github.com/prabhakarlab/Banksy/)

``` r
library(Banksy)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)

library(ggplot2)
library(gridExtra)
library(pals)

# Kelly palette for visualization
mypal <- kelly()[-1]
```

## Running BANKSY within Seurat’s spatial framework

We demonstrate how to run BANKSY within Seurat’s spatial analysis
framework with a mouse hippocampus Slide-seq v2 dataset from the
*SeuratData* package.

After installing *SeuratData*, the data can be accessed as follows:

``` r
InstallData('ssHippo')
ss.hippo <- LoadData("ssHippo")
```

We perform simple preprocessing by filtering beads with high mito
percentage and keeping only beads within the 5th and 98th percentile of
total UMI counts. To keep runtime of this vignette short, we downsample
the data to 10,000 beads.

``` r
# Filtering
ss.hippo[['percent.mt']] <- PercentageFeatureSet(ss.hippo, pattern = '^MT-')
ss.hippo <- subset(ss.hippo, percent.mt < 10 &
                    nCount_Spatial > quantile(ss.hippo$nCount_Spatial, 0.05) &
                    nCount_Spatial < quantile(ss.hippo$nCount_Spatial, 0.98))
# Downsample
set.seed(42)
ss.hippo <- ss.hippo[,sample(colnames(ss.hippo), 1e4)]
```

Next, normalize the data and find variable genes:

``` r
# Normalize
ss.hippo <- NormalizeData(ss.hippo)
ss.hippo <- FindVariableFeatures(ss.hippo)
ss.hippo <- ScaleData(ss.hippo)
```

To run BANKSY, we specify the following:

- `lambda`: a numeric value in \[0,1\]. With low values of lambda,
  BANKSY operates in cell-typing mode, while high values of lambda find
  spatial domains.
- `assay` and `slot`: determines where to pull the expression data from
- `features`: specifies features for downstream analysis. This can be
  `'all'`, `'variable'` or a subset of features.  
- `k_geom`: the number of neighbors that defines a cell’s neighborhood

Call `?RunBanksy` for more details on function parameters.

``` r
# Run BANKSY
ss.hippo <- RunBanksy(ss.hippo, lambda = 0.2, verbose=TRUE, 
                      assay = 'Spatial', slot = 'data', features = 'variable',
                      k_geom = 15)
ss.hippo
```

    ## An object of class Seurat 
    ## 27264 features across 10000 samples within 2 assays 
    ## Active assay: BANKSY (4000 features, 0 variable features)
    ##  2 layers present: data, scale.data
    ##  1 other assay present: Spatial
    ##  1 image present: image

Note that the `RunBanksy` function sets the default assay to `BANKSY` (
determined by the `assay_name` argument) and fills the `scale.data`
slot. Users should not call `ScaleData` on the `BANKSY` assay as this
negates the effects of `lambda`.

The rest of the pipeline is similar to the ‘default’ Seurat pipline. We
scale the data and run dimensionality reduction with PCA and UMAP:

``` r
# Run PCA and UMAP
ss.hippo <- RunPCA(ss.hippo, assay = 'BANKSY', features = rownames(ss.hippo), npcs = 30)
ss.hippo <- RunUMAP(ss.hippo, dims = 1:30)
```

Next, find BANKSY clusters:

``` r
# Clustering
ss.hippo <- FindNeighbors(ss.hippo, dims = 1:30)
ss.hippo <- FindClusters(ss.hippo, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10000
    ## Number of edges: 365658
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9033
    ## Number of communities: 13
    ## Elapsed time: 1 seconds

Visualize the UMAP and Spatial plot:

``` r
# Viz
grid.arrange(
    DimPlot(ss.hippo, pt.size = 0.25, label = TRUE, label.size = 3, repel = TRUE),
    SpatialDimPlot(ss.hippo, stroke = NA, label = TRUE, label.size = 3, 
                   repel = TRUE, alpha = 0.5, pt.size.factor = 2),
    ncol = 2
)
```

<img src="banksy_files/figure-gfm/ss_viz-1.png" style="display: block; margin: auto;" />

Find markers based on the BANKSY clusters and visualize them. Here, we
find differentially expressed genes between the CA1 and CA3 regions.

``` r
# Find markers
DefaultAssay(ss.hippo) <- 'Spatial'
markers <- FindMarkers(ss.hippo, ident.1 = 4, ident.2 = 9, only.pos = F, 
                       logfc.threshold = 1, min.pct = 0.5)
markers <- markers[markers$p_val_adj < 0.01,]
markers
```

    ##               p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## SNAP25 1.127235e-46  -1.260312 0.658 0.823 2.622400e-42
    ## CHGB   9.840001e-44  -1.985343 0.439 0.697 2.289178e-39
    ## STMN2  1.281230e-24  -1.430138 0.335 0.574 2.980653e-20
    ## SYN2   3.272800e-23  -1.609355 0.332 0.564 7.613842e-19
    ## ATP2B1 1.545647e-22   1.251540 0.639 0.474 3.595793e-18
    ## CPLX2  4.619232e-21  -1.220110 0.289 0.522 1.074618e-16
    ## PRKCB  1.276453e-18   1.394809 0.552 0.341 2.969539e-14
    ## PCP4   2.006224e-18  -1.269671 0.379 0.578 4.667279e-14
    ## TUBB2A 1.330787e-16  -1.054176 0.450 0.629 3.095942e-12
    ## DDN    1.784378e-14   1.401976 0.592 0.396 4.151176e-10
    ## SNCA   7.596526e-12  -1.022314 0.397 0.544 1.767256e-07

``` r
genes <- c('ATP2B1', 'CHGB')
SpatialFeaturePlot(ss.hippo, features = genes, pt.size.factor = 3, 
                   stroke = NA, alpha = 0.5, max.cutoff = 'q95')
```

<img src="banksy_files/figure-gfm/ss_markers-1.png" style="display: block; margin: auto;" />

## Running BANKSY with locations provided explicitly

One can also call `RunBanksy` on a Seurat object created from counts by
providing the location of cell centroids or spots explicitly. In this
case, the locations must be stored as metadata. Here, we use a mouse
hippocampus VeraFISH dataset provided with the *Banksy* package.

``` r
data(hippocampus)
head(hippocampus$expression[,1:5])
```

    ##         cell_1276 cell_8890 cell_691 cell_396 cell_9818
    ## Sparcl1        45         0       11       22         0
    ## Slc1a2         17         0        6        5         0
    ## Map            10         0       12       16         0
    ## Sqstm1         26         0        0        2         0
    ## Atp1a2          0         0        4        3         0
    ## Tnc             0         0        0        0         0

``` r
head(hippocampus$locations)
```

    ##                 sdimx    sdimy
    ## cell_1276  -13372.899 15776.37
    ## cell_8890    8941.101 15866.37
    ## cell_691   -14882.899 15896.37
    ## cell_396   -15492.899 15835.37
    ## cell_9818   11308.101 15846.37
    ## cell_11310  14894.101 15810.37

Construct the Seurat object by storing the locations of cell centroids
as metadata. We keep cells with total count between 5th and 98th
percentile:

``` r
# Create manually
vf.hippo <- CreateSeuratObject(counts = hippocampus$expression,
                               meta.data = hippocampus$locations)
vf.hippo <- subset(vf.hippo,
                   nCount_RNA > quantile(vf.hippo$nCount_RNA, 0.05) & 
                   nCount_RNA < quantile(vf.hippo$nCount_RNA, 0.98))
```

Next, we normalize the data by library size and scale the data:

``` r
# Normalize
vf.hippo <- NormalizeData(vf.hippo, scale.factor = 100, normalization.method = 'RC')
vf.hippo <- ScaleData(vf.hippo)
```

Now, run BANKSY. Here, we provide the column names of the x and y
spatial coordinates as stored in the metadata to `dimx` and `dimy`
respectively:

``` r
# Run BANKSY
vf.hippo <- RunBanksy(vf.hippo, lambda = 0.2, dimx = 'sdimx', dimy = 'sdimy', 
                      assay = 'RNA', slot = 'data', features = 'all', k_geom = 10)
```

Note that the `RunBanksy` function sets the default assay to `BANKSY` (
determined by the `assay_name` argument) and fills the `scale.data`
slot. Users should not call `ScaleData` on the `BANKSY` assay as this
negates the effects of `lambda`.

Run PCA on the BANKSY matrix:

``` r
# PCA
vf.hippo <- RunPCA(vf.hippo, assay = 'BANKSY', features = rownames(vf.hippo), npcs = 20)
```

Find BANKSY clusters:

``` r
# Cluster
vf.hippo <- FindNeighbors(vf.hippo, dims = 1:20)
vf.hippo <- FindClusters(vf.hippo, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10205
    ## Number of edges: 446178
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9099
    ## Number of communities: 15
    ## Elapsed time: 1 seconds

Visualise BANKSY clusters in spatial dimensions:

``` r
# Viz
FeatureScatter(vf.hippo, 'sdimx', 'sdimy', cols = mypal, pt.size = 0.75)
```

<img src="banksy_files/figure-gfm/hippo_viz-1.png" style="display: block; margin: auto;" />

``` r
FeatureScatter(vf.hippo, 'sdimx', 'sdimy', cols = mypal, pt.size = 0.1) + facet_wrap(~ colors)
```

<img src="banksy_files/figure-gfm/hippo_viz-2.png" style="display: block; margin: auto;" />

Find markers and visualise them. Here, we do so for a cluster defined by
a thin layer of cells expressing Gfap. We also write a simple function
`genePlot` that plots marker genes in spatial dimensions.

``` r
# Find markers
DefaultAssay(vf.hippo) <- 'RNA'
markers <- FindMarkers(vf.hippo, ident.1 = 6, only.pos = TRUE)

genePlot <- function(object, dimx, dimy, gene, assay = 'RNA',
                     slot = 'scale.data', q.low = 0.01, q.high = 0.99,
                     col.low='blue', col.high='red') {
    val <- GetAssayData(object, assay=assay, slot=slot)[gene,]
    val.low <- quantile(val, q.low)
    val.high <- quantile(val, q.high)
    val[val < val.low] <- val.low
    val[val > val.high] <- val.high
    pdf <- data.frame(x=object[[dimx]], y=object[[dimy]], gene=val)
    colnames(pdf) <- c('sdimx','sdimy', 'gene')
    ggplot(pdf, aes(x=sdimx,y=sdimy,color=gene)) + geom_point(size = 1) + 
        theme_minimal() + theme(legend.title = element_blank()) +
        scale_color_gradient2(low = col.low, high = col.high) +
        ggtitle(gene)
}

genePlot(vf.hippo, 'sdimx', 'sdimy', 'Gfap')
```

<img src="banksy_files/figure-gfm/hippo_gene-1.png" style="display: block; margin: auto;" />

## Multi-sample analysis

This section demonstrate demonstrates multi-sample analysis. Such an
approach is appropriate when analysing multiple spatial omics datasets
with non-contiguous spatial coordinates, and when large batch effects
are not present.

Here, we use a mouse hippocampus VeraFISH dataset provided with the
*Banksy* package.

``` r
data(hippocampus)
head(hippocampus$expression[,1:5])
```

    ##         cell_1276 cell_8890 cell_691 cell_396 cell_9818
    ## Sparcl1        45         0       11       22         0
    ## Slc1a2         17         0        6        5         0
    ## Map            10         0       12       16         0
    ## Sqstm1         26         0        0        2         0
    ## Atp1a2          0         0        4        3         0
    ## Tnc             0         0        0        0         0

``` r
head(hippocampus$locations)
```

    ##                 sdimx    sdimy
    ## cell_1276  -13372.899 15776.37
    ## cell_8890    8941.101 15866.37
    ## cell_691   -14882.899 15896.37
    ## cell_396   -15492.899 15835.37
    ## cell_9818   11308.101 15846.37
    ## cell_11310  14894.101 15810.37

For demonstration purposes, we create three separate datasets by
splitting the data.

``` r
# Number of groups
n_groups = 3
group_names = paste0('group', seq(n_groups))
group_size = 1000
starts = seq(1, by=group_size, length.out=n_groups)
ends = starts + group_size - 1

# List of Seurat objects
seu_list = lapply(seq(n_groups), function(i) {
  idx = seq(starts[i], ends[i])
  seu = CreateSeuratObject(
    counts = hippocampus$expression[,idx],
    meta.data = data.frame(scale(hippocampus$locations[idx,], scale = FALSE))
  )
  # Set original identity of cell
  seu$orig.ident = group_names[i]
  seu
})
seu_list
```

    ## [[1]]
    ## An object of class Seurat 
    ## 120 features across 1000 samples within 1 assay 
    ## Active assay: RNA (120 features, 0 variable features)
    ##  1 layer present: counts
    ## 
    ## [[2]]
    ## An object of class Seurat 
    ## 120 features across 1000 samples within 1 assay 
    ## Active assay: RNA (120 features, 0 variable features)
    ##  1 layer present: counts
    ## 
    ## [[3]]
    ## An object of class Seurat 
    ## 120 features across 1000 samples within 1 assay 
    ## Active assay: RNA (120 features, 0 variable features)
    ##  1 layer present: counts

Perform normalisation for each dataset.

``` r
seu_list = lapply(seu_list, NormalizeData,
                  scale.factor = 100, normalization.method = 'RC')
```

Merge the datasets. Note that the spatial coordinates overlap.

``` r
# Merge
seu = Reduce(merge, seu_list)
seu = JoinLayers(seu) # run this for Seurat v5 objects

# Plot spatial coordinates colored by group
plot(FetchData(seu, c('sdimx', 'sdimy')), col = factor(seu$orig.ident))
```

<img src="banksy_files/figure-gfm/multi-spatial-1.png" style="display: block; margin: auto;" />

Now run BANKSY. For multi-sample analysis, the argument `group` must be
provided, which specifies the name of the metadata column that gives the
assignment of each cell or spot to its original Seurat object. Here, we
use `orig.ident`. Internally, providing the `group` argument tells the
function to compute neighborhood matrices based on locations staggered
by `group`, ensuring that cells from different spatial datasets do not
overlap. The staggered locations are stored in the metadata for sanity
checking. The `split.scale` argument allows for within-group scaling,
accounting for minor differences in datasets.

``` r
# Grouping variable
head(seu@meta.data)
```

    ##            orig.ident nCount_RNA nFeature_RNA     sdimx    sdimy
    ## cell_1276      group1        266           51 -11933.19 1366.934
    ## cell_8890      group1         13            3  10380.81 1456.934
    ## cell_691       group1        132           36 -13443.19 1486.934
    ## cell_396       group1         95           27 -14053.19 1425.934
    ## cell_9818      group1         10            5  12747.81 1436.934
    ## cell_11310     group1         15            5  16333.81 1400.934

``` r
table(seu$orig.ident)
```

    ## 
    ## group1 group2 group3 
    ##   1000   1000   1000

``` r
# Run BANKSY
seu = RunBanksy(seu, lambda = 0.2, assay = 'RNA', slot = 'data',
                dimx = 'sdimx', dimy = 'sdimy', features = 'all',
                group = 'orig.ident', split.scale = TRUE, k_geom = 15)

# Staggered locations added to metadata
head(seu@meta.data)
```

    ##            orig.ident nCount_RNA nFeature_RNA     sdimx    sdimy
    ## cell_1276      group1        266           51 -11933.19 1366.934
    ## cell_8890      group1         13            3  10380.81 1456.934
    ## cell_691       group1        132           36 -13443.19 1486.934
    ## cell_396       group1         95           27 -14053.19 1425.934
    ## cell_9818      group1         10            5  12747.81 1436.934
    ## cell_11310     group1         15            5  16333.81 1400.934
    ##            staggered_sdimx staggered_sdimy
    ## cell_1276         3728.686        1366.934
    ## cell_8890        26042.686        1456.934
    ## cell_691          2218.686        1486.934
    ## cell_396          1608.686        1425.934
    ## cell_9818        28409.686        1436.934
    ## cell_11310       31995.686        1400.934

The rest of the workflow follows as before:

``` r
seu = RunPCA(seu, assay = 'BANKSY', features = rownames(seu), npcs = 30)
seu = RunUMAP(seu, dims = 1:30)
seu = FindNeighbors(seu, dims = 1:30)
seu = FindClusters(seu, resolution = 1)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3000
    ## Number of edges: 171757
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8094
    ## Number of communities: 12
    ## Elapsed time: 0 seconds

Visualise clusters:

``` r
mypal <- kelly()[-1]
DimPlot(seu, pt.size = 0.25, label = TRUE, label.size = 3, cols = mypal)
```

<img src="banksy_files/figure-gfm/multi-umap-1.png" style="display: block; margin: auto;" />

``` r
FeatureScatter(seu, 'staggered_sdimx', 'staggered_sdimy', pt.size = 0.75, cols = mypal)
```

<img src="banksy_files/figure-gfm/multi-spatial-staggered-1.png" style="display: block; margin: auto;" />

## Spatial data integration with Harmony

BANKSY can be used with Harmony for integrating multiple spatial omics
datasets in the presence of strong batch effects.

Download the data.

``` r
library(spatialLIBD)
library(ExperimentHub)
library(harmony)

ehub <- ExperimentHub::ExperimentHub()
spe <- spatialLIBD::fetch_data(type = "spe", eh = ehub)

imgData(spe) <- NULL
assay(spe, "logcounts") <- NULL
reducedDims(spe) <- NULL
rowData(spe) <- NULL
colData(spe) <- DataFrame(
  sample_id = spe$sample_id,
  clust_annotation = factor(
    addNA(spe$layer_guess_reordered_short),
    exclude = NULL, labels = seq(8)
  ),
  in_tissue = spe$in_tissue,
  row.names = colnames(spe)
)
invisible(gc())

# Subset to first sample of each subject
sample_names <- c("151507", "151669", "151673")
spe_list <- lapply(sample_names, function(x) spe[, spe$sample_id == x])
rm(spe)
invisible(gc())
```

Normalise the data and compute highly variable features.

``` r
# Convert to Seurat and Normalize data
seu_list <- lapply(spe_list, function(x) {
  x <- as.Seurat(x, data = NULL)
  NormalizeData(x, scale.factor = 3000, normalization.method = 'RC')
})

# Compute HVGs for each dataset and take the union
hvgs <- lapply(seu_list, function(x) {
  VariableFeatures(FindVariableFeatures(x, nfeatures = 2000))
})
hvgs <- Reduce(union, hvgs)

# Subset to HVGs
seu_list <- lapply(seu_list, function(x) x[hvgs,])
seu <- Reduce(merge, seu_list)

locs <- do.call(rbind.data.frame, lapply(spe_list, spatialCoords))
seu@meta.data <- cbind(seu@meta.data, locs)

seu
```

Run BANKSY. When analysing multiple samples, the argument `group` must
be provided, which specifies the name of the metadata column that gives
the assignment of each cell or spot to its original Seurat object. Here,
we use `sample_id`. Internally, providing the `group` argument tells the
function to compute neighborhood matrices based on locations staggered
by `group`, ensuring that cells from different spatial datasets do not
overlap. The staggered locations are stored in the metadata for sanity
checking. Within-group scaling has little effect in the presence of
strong batch effects, hence, we set `split.scale=FALSE` for efficiency.

``` r
# Grouping variable
head(seu@meta.data)
table(seu$sample_id)

sdimx <- 'pxl_col_in_fullres'
sdimy <- 'pxl_row_in_fullres'

# Run BANKSY
seu <- RunBanksy(seu, lambda = 0.2, assay = 'originalexp', slot = 'data',
                dimx = sdimx, dimy = sdimy, features = 'all',
                group = 'sample_id', split.scale = FALSE, k_geom = 6)
```

Compute a spatially-aware embedding with PCA on the BANKSY matrix, and
run Harmony on this embedding.

``` r
seu <- RunPCA(seu, assay = 'BANKSY', features = rownames(seu), npcs = 10)
seu <- RunHarmony(seu, group.by.vars='sample_id')
```

The rest of the workflow follows as before:

``` r
seu <- RunUMAP(seu, dims = 1:10, reduction = 'harmony')
seu <- FindNeighbors(seu, dims = 1:10, reduction = 'harmony')
seu <- FindClusters(seu, resolution = 0.4)
```

Visualise clusters:

``` r
DimPlot(seu, pt.size = 0.25, label = TRUE, label.size = 3, cols = mypal)
FeatureScatter(seu, 'staggered_sdimx', 'staggered_sdimy', cols = mypal, pt.size = 0.75)
```

## Getting help

For more information, visit <https://github.com/prabhakarlab/Banksy>.

<details>
<summary>
Vignette runtime
</summary>

    ## Time difference of 1.434785 mins

</details>
<details>
<summary>
Session info
</summary>

``` r
sessionInfo()
```

    ## R version 4.3.2 (2023-10-31)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Sonoma 14.2.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pals_1.8                 gridExtra_2.3            ggplot2_3.4.4           
    ##  [4] SeuratWrappers_0.3.4     ssHippo.SeuratData_3.1.4 SeuratData_0.2.2.9001   
    ##  [7] Seurat_5.0.1             SeuratObject_5.0.1       sp_2.1-3                
    ## [10] Banksy_0.99.9           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppHungarian_0.3           RcppAnnoy_0.0.22           
    ##   [3] splines_4.3.2               later_1.3.2                
    ##   [5] bitops_1.0-7                tibble_3.2.1               
    ##   [7] R.oo_1.26.0                 polyclip_1.10-6            
    ##   [9] fastDummies_1.7.3           lifecycle_1.0.4            
    ##  [11] aricode_1.0.3               globals_0.16.2             
    ##  [13] lattice_0.22-5              MASS_7.3-60.0.1            
    ##  [15] magrittr_2.0.3              limma_3.58.1               
    ##  [17] plotly_4.10.4               rmarkdown_2.25             
    ##  [19] yaml_2.3.8                  remotes_2.4.2.1            
    ##  [21] httpuv_1.6.14               sctransform_0.4.1          
    ##  [23] spam_2.10-0                 spatstat.sparse_3.0-3      
    ##  [25] reticulate_1.35.0           mapproj_1.2.11             
    ##  [27] cowplot_1.1.3               pbapply_1.7-2              
    ##  [29] RColorBrewer_1.1-3          maps_3.4.2                 
    ##  [31] abind_1.4-5                 zlibbioc_1.48.0            
    ##  [33] Rtsne_0.17                  GenomicRanges_1.54.1       
    ##  [35] purrr_1.0.2                 R.utils_2.12.3             
    ##  [37] BiocGenerics_0.48.1         RCurl_1.98-1.14            
    ##  [39] rappdirs_0.3.3              GenomeInfoDbData_1.2.11    
    ##  [41] IRanges_2.36.0              S4Vectors_0.40.2           
    ##  [43] ggrepel_0.9.5               irlba_2.3.5.1              
    ##  [45] listenv_0.9.1               spatstat.utils_3.0-4       
    ##  [47] goftest_1.2-3               RSpectra_0.16-1            
    ##  [49] spatstat.random_3.2-2       fitdistrplus_1.1-11        
    ##  [51] parallelly_1.37.0           leiden_0.4.3.1             
    ##  [53] codetools_0.2-19            DelayedArray_0.28.0        
    ##  [55] tidyselect_1.2.0            farver_2.1.1               
    ##  [57] matrixStats_1.2.0           stats4_4.3.2               
    ##  [59] spatstat.explore_3.2-6      jsonlite_1.8.8             
    ##  [61] ellipsis_0.3.2              progressr_0.14.0           
    ##  [63] ggridges_0.5.6              survival_3.5-7             
    ##  [65] dbscan_1.1-12               tools_4.3.2                
    ##  [67] ica_1.0-3                   Rcpp_1.0.12                
    ##  [69] glue_1.7.0                  SparseArray_1.2.4          
    ##  [71] xfun_0.42                   MatrixGenerics_1.14.0      
    ##  [73] GenomeInfoDb_1.38.6         dplyr_1.1.4                
    ##  [75] withr_3.0.0                 BiocManager_1.30.22        
    ##  [77] fastmap_1.1.1               fansi_1.0.6                
    ##  [79] digest_0.6.34               rsvd_1.0.5                 
    ##  [81] R6_2.5.1                    mime_0.12                  
    ##  [83] colorspace_2.1-0            scattermore_1.2            
    ##  [85] sccore_1.0.4                tensor_1.5                 
    ##  [87] dichromat_2.0-0.1           spatstat.data_3.0-4        
    ##  [89] R.methodsS3_1.8.2           utf8_1.2.4                 
    ##  [91] tidyr_1.3.1                 generics_0.1.3             
    ##  [93] data.table_1.15.0           httr_1.4.7                 
    ##  [95] htmlwidgets_1.6.4           S4Arrays_1.2.0             
    ##  [97] uwot_0.1.16                 pkgconfig_2.0.3            
    ##  [99] gtable_0.3.4                lmtest_0.9-40              
    ## [101] SingleCellExperiment_1.24.0 XVector_0.42.0             
    ## [103] htmltools_0.5.7             dotCall64_1.1-1            
    ## [105] scales_1.3.0                Biobase_2.62.0             
    ## [107] png_0.1-8                   SpatialExperiment_1.12.0   
    ## [109] knitr_1.45                  rstudioapi_0.15.0          
    ## [111] reshape2_1.4.4              rjson_0.2.21               
    ## [113] nlme_3.1-164                zoo_1.8-12                 
    ## [115] stringr_1.5.1               KernSmooth_2.23-22         
    ## [117] parallel_4.3.2              miniUI_0.1.1.1             
    ## [119] pillar_1.9.0                grid_4.3.2                 
    ## [121] vctrs_0.6.5                 RANN_2.6.1                 
    ## [123] promises_1.2.1              xtable_1.8-4               
    ## [125] cluster_2.1.6               evaluate_0.23              
    ## [127] magick_2.8.2                cli_3.6.2                  
    ## [129] compiler_4.3.2              rlang_1.1.3                
    ## [131] crayon_1.5.2                future.apply_1.11.1        
    ## [133] labeling_0.4.3              mclust_6.0.1               
    ## [135] plyr_1.8.9                  stringi_1.8.3              
    ## [137] viridisLite_0.4.2           deldir_2.0-2               
    ## [139] munsell_0.5.0               lazyeval_0.2.2             
    ## [141] spatstat.geom_3.2-8         Matrix_1.6-5               
    ## [143] RcppHNSW_0.6.0              patchwork_1.2.0            
    ## [145] future_1.33.1               statmod_1.5.0              
    ## [147] shiny_1.8.0                 highr_0.10                 
    ## [149] SummarizedExperiment_1.32.0 ROCR_1.0-11                
    ## [151] leidenAlg_1.1.2             igraph_2.0.1.1

</details>
