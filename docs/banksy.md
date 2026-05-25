Running BANKSY with Seurat
================
Compiled: May 25, 2026

- [Introduction](#introduction)
- [Overview](#overview)
- [Running BANKSY within Seurat’s spatial
  framework](#running-banksy-within-seurats-spatial-framework)
- [Running BANKSY with locations provided
  explicitly](#running-banksy-with-locations-provided-explicitly)
- [Multi-sample analysis](#multi-sample-analysis)
- [Spatial data integration with
  Harmony](#spatial-data-integration-with-harmony)
- [Scaling to large datasets](#scaling-to-large-datasets)
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
  `'all'`, `'variable'` or a subset of features.\
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
    ## Number of edges: 366223
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9036
    ## Number of communities: 13
    ## Elapsed time: 0 seconds

Visualize the UMAP and Spatial plot:

``` r
# Viz
coords <- GetTissueCoordinates(ss.hippo)
coords$cluster <- Idents(ss.hippo)[coords$cells]
spatial_plot <- ggplot(coords, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 0.25) + theme_minimal() +
    guides(color = guide_legend(override.aes = list(size = 3)))
grid.arrange(
    DimPlot(ss.hippo, pt.size = 0.25, label = TRUE, label.size = 3, repel = TRUE),
    spatial_plot,
    ncol = 2
)
```

<img src="banksy_files/figure-gfm/ss_viz-1.png" alt="" style="display: block; margin: auto;" />

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
    ## SNAP25 5.201207e-47  -1.259593 0.657 0.825 1.210009e-42
    ## CHGB   2.606034e-44  -1.991969 0.436 0.699 6.062677e-40
    ## STMN2  2.534375e-25  -1.448304 0.333 0.576 5.895969e-21
    ## SYN2   6.294296e-23  -1.566838 0.333 0.564 1.464305e-18
    ## ATP2B1 2.152782e-22   1.249215 0.638 0.474 5.008231e-18
    ## CPLX2  9.743390e-22  -1.240850 0.287 0.524 2.266702e-17
    ## PRKCB  1.614673e-18   1.392484 0.551 0.341 3.756375e-14
    ## PCP4   1.680405e-18  -1.271996 0.378 0.578 3.909294e-14
    ## TUBB2A 1.107377e-16  -1.056501 0.449 0.629 2.576203e-12
    ## DDN    2.213579e-14   1.399651 0.591 0.396 5.149670e-10
    ## SNCA   3.187591e-12  -1.031572 0.395 0.546 7.415612e-08

``` r
genes <- c('ATP2B1', 'CHGB')
coords <- GetTissueCoordinates(ss.hippo)
gene_plots <- lapply(genes, function(g) {
    val <- GetAssayData(ss.hippo, assay = 'Spatial', layer = 'data')[g, coords$cells]
    val[val > quantile(val, 0.95)] <- quantile(val, 0.95)
    coords$expr <- val
    ggplot(coords, aes(x = x, y = y, color = expr)) +
        geom_point(size = 0.25) + theme_minimal() +
        scale_color_gradient(low = 'grey90', high = 'red') +
        ggtitle(g) + theme(legend.title = element_blank())
})
grid.arrange(grobs = gene_plots, ncol = 2)
```

<img src="banksy_files/figure-gfm/ss_markers-1.png" alt="" style="display: block; margin: auto;" />

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
    ## Elapsed time: 0 seconds

Visualise BANKSY clusters in spatial dimensions:

``` r
# Viz
FeatureScatter(vf.hippo, 'sdimx', 'sdimy', cols = mypal, pt.size = 0.75)
```

<img src="banksy_files/figure-gfm/hippo_viz-1.png" alt="" style="display: block; margin: auto;" />

``` r
FeatureScatter(vf.hippo, 'sdimx', 'sdimy', cols = mypal, pt.size = 0.1) + facet_wrap(~ colors)
```

<img src="banksy_files/figure-gfm/hippo_viz-2.png" alt="" style="display: block; margin: auto;" />

Find markers and visualise them. Here, we do so for a cluster defined by
a thin layer of cells expressing Gfap. We also write a simple function
`genePlot` that plots marker genes in spatial dimensions.

``` r
# Find markers
DefaultAssay(vf.hippo) <- 'RNA'
markers <- FindMarkers(vf.hippo, ident.1 = 6, only.pos = TRUE)

genePlot <- function(object, dimx, dimy, gene, assay = 'RNA',
                     layer = 'scale.data', q.low = 0.01, q.high = 0.99,
                     col.low='blue', col.high='red') {
    val <- GetAssayData(object, assay=assay, layer=layer)[gene,]
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

<img src="banksy_files/figure-gfm/hippo_gene-1.png" alt="" style="display: block; margin: auto;" />

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

<img src="banksy_files/figure-gfm/multi-spatial-1.png" alt="" style="display: block; margin: auto;" />

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

<img src="banksy_files/figure-gfm/multi-umap-1.png" alt="" style="display: block; margin: auto;" />

``` r
FeatureScatter(seu, 'staggered_sdimx', 'staggered_sdimy', pt.size = 0.75, cols = mypal)
```

<img src="banksy_files/figure-gfm/multi-spatial-staggered-1.png" alt="" style="display: block; margin: auto;" />

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

    ## An object of class Seurat 
    ## 4467 features across 11526 samples within 1 assay 
    ## Active assay: originalexp (4467 features, 0 variable features)
    ##  2 layers present: counts, data

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
```

    ##                         orig.ident nCount_originalexp nFeature_originalexp
    ## AAACAACGAATAGTTC-1_1 SeuratProject                288                  197
    ## AAACAAGTATCTCCCA-1_1 SeuratProject               1546                  550
    ## AAACAATCTACTAGCA-1_1 SeuratProject                525                  247
    ## AAACACCAATAACTGC-1_1 SeuratProject               1123                  457
    ## AAACAGCTTTCAGAAG-1_1 SeuratProject               1008                  415
    ## AAACAGGGTCTATATT-1_1 SeuratProject               1431                  539
    ##                      sample_id clust_annotation in_tissue pxl_col_in_fullres
    ## AAACAACGAATAGTTC-1_1    151507                1      TRUE               3276
    ## AAACAAGTATCTCCCA-1_1    151507                3      TRUE               9178
    ## AAACAATCTACTAGCA-1_1    151507                1      TRUE               5133
    ## AAACACCAATAACTGC-1_1    151507                7      TRUE               3462
    ## AAACAGCTTTCAGAAG-1_1    151507                6      TRUE               2779
    ## AAACAGGGTCTATATT-1_1    151507                6      TRUE               3053
    ##                      pxl_row_in_fullres
    ## AAACAACGAATAGTTC-1_1               2514
    ## AAACAAGTATCTCCCA-1_1               8520
    ## AAACAATCTACTAGCA-1_1               2878
    ## AAACACCAATAACTGC-1_1               9581
    ## AAACAGCTTTCAGAAG-1_1               7663
    ## AAACAGGGTCTATATT-1_1               8143

``` r
table(seu$sample_id)
```

    ## 
    ## 151507 151669 151673 
    ##   4226   3661   3639

``` r
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

    ## Initializing centroids

The rest of the workflow follows as before:

``` r
seu <- RunUMAP(seu, dims = 1:10, reduction = 'harmony')
seu <- FindNeighbors(seu, dims = 1:10, reduction = 'harmony')
seu <- FindClusters(seu, resolution = 0.4)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 11526
    ## Number of edges: 373678
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8790
    ## Number of communities: 7
    ## Elapsed time: 1 seconds

Visualise clusters:

``` r
grid.arrange(
    DimPlot(seu, pt.size = 0.25, label = TRUE, label.size = 3, cols = mypal),
    DimPlot(seu, pt.size = 0.25, group.by = 'sample_id'),
    ncol = 2
)
```

<img src="banksy_files/figure-gfm/harmony_viz-1.png" alt="" style="display: block; margin: auto;" />

``` r
spatial_plots <- lapply(sample_names, function(s) {
    sub <- subset(seu, sample_id == s)
    ggplot(sub@meta.data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
                              color = seurat_clusters)) +
        geom_point(size = 1) + coord_equal() + theme_minimal() +
        scale_color_manual(values = mypal) +
        ggtitle(s) + theme(legend.position = 'none') +
        labs(x = NULL, y = NULL)
})
grid.arrange(grobs = spatial_plots, ncol = 3)
```

<img src="banksy_files/figure-gfm/harmony_spatial-1.png" alt="" style="display: block; margin: auto;" />

## Scaling to large datasets

For large datasets (millions of cells), materializing the full BANKSY
matrix can exceed memory limits. The `lazy` option in `RunBanksy` avoids
this by computing PCs of the BANKSY matrix directly via a lazy linear
operator. The BANKSY matrix is never formed; instead, matrix-vector
products are evaluated on-the-fly using the sparse expression matrix and
weight matrix. This reduces peak memory from
$O(\text{genes} \times \text{cells})$ (the dense BANKSY matrix) to
$O(\text{nnz}(\text{input}))$ (the number of non-zero entries in the
sparse input), since only the sparse input and sparse weight matrix are
held in memory.

**Note**: `lazy=TRUE` currently supports M=0 only (no azimuthal Gabor
filter). When `group` is specified, `split.scale=TRUE` applies
within-group scaling in the lazy workflow as in the standard workflow.

This feature requires the development versions of Banksy and
SeuratWrappers:

``` r
devtools::install_github('prabhakarlab/Banksy@feat-sparse-matmul')
devtools::install_github('jleechung/seurat-wrappers@feat-sparse-matmul')
```

We demonstrate this with a 10X Xenium mouse brain dataset (~37k cells),
performing domain segmentation with `lambda = 0.8`.

Download the dataset:

``` r
download.file(
    url = 'https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip',
    destfile = '~/Downloads/xenium_subset.zip'
)
unzip('~/Downloads/xenium_subset.zip', exdir = '~/Downloads/xenium_subset')
```

Load and preprocess:

``` r
xenium <- LoadXenium('~/Downloads/xenium_subset', fov = 'fov')
xenium <- NormalizeData(xenium)
xenium <- FindVariableFeatures(xenium)
```

With `lazy = TRUE`, `RunBanksy` reads the expression matrix from the
`data` slot of the `Xenium` assay, computes `npcs` principal components
of the BANKSY matrix (default 50, matching Seurat’s `RunPCA`) directly
via a sparse linear operator, and stores the result as a dimensionality
reduction named by `assay_name`. No `BANKSY` assay is created, and there
is no need to call `RunPCA` separately. For multi-sample integration,
you can run Harmony on the BANKSY reduction before `FindNeighbors` (see
the [Harmony section](#spatial-data-integration-with-harmony) above).

``` r
xenium_input <- xenium
xenium <- RunBanksy(xenium, lazy = TRUE, npcs = 50,
                    lambda = 0.8, k_geom = 30, use_agf = FALSE,
                    assay = 'Xenium', slot = 'data', features = 'variable',
                    assay_name = 'BANKSY', verbose = TRUE)

# Optionally, run Harmony for batch correction
# xenium <- RunHarmony(xenium, group.by.vars = 'batch', reduction.use = 'BANKSY')

xenium <- FindNeighbors(xenium, reduction = 'BANKSY', dims = 1:50)
xenium <- FindClusters(xenium, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 36602
    ## Number of edges: 1017035
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9441
    ## Number of communities: 23
    ## Elapsed time: 3 seconds

Visualize the spatial domains:

``` r
n_clust <- length(levels(Idents(xenium)))
pal <- c(mypal, hcl.colors(max(0, n_clust - length(mypal))))
ImageDimPlot(xenium, size = 0.5, cols = pal) +
    coord_flip() + scale_x_reverse() + NoLegend() +
    theme(plot.background = element_rect(fill = 'white'),
          panel.background = element_rect(fill = 'white'))
```

<img src="banksy_files/figure-gfm/xenium_viz-1.png" alt="" style="display: block; margin: auto;" />

<details>

<summary>

Equivalence with the standard workflow
</summary>

The standard workflow constructs the full BANKSY matrix `M` by
vertically stacking the z-scored expression matrix and the z-scored
neighborhood-smoothed matrix (weighted by `lambda`), where z-scoring
includes clipping to \[-10, 10\] to match Seurat’s `FastRowScale`. PCA
is then computed on `M` via `RunPCA`.

With `lazy=TRUE`, `M` is never formed. Instead, `irlba` accesses `M`
through a lazy linear operator that evaluates matrix-vector products on
the fly. The key insight is that `(XW)v = X(Wv)` by associativity: the
sparse neighbor-weight matrix `W` is applied to the vector first, then
left-multiplied by the sparse expression matrix `X`, avoiding formation
of the dense product `XW`. Row means, standard deviations, and a sparse
correction for z-score clipping are precomputed once before the Lanczos
iterations begin.

This reduces peak memory from `O(g*n)` (the dense BANKSY matrix) to
`O(nnz(X) + kn)` (the sparse input plus the sparse weight matrix), while
producing numerically identical PCA embeddings.

The standard workflow materializes the full BANKSY matrix as a Seurat
assay, then runs PCA via `RunPCA`:

``` r
xenium_std <- RunBanksy(xenium_input, lambda = 0.8, verbose = TRUE,
                        assay = 'Xenium', slot = 'data', features = 'variable',
                        use_agf = FALSE, k_geom = 30)
xenium_std <- RunPCA(xenium_std, assay = 'BANKSY', features = rownames(xenium_std),
                     npcs = 50)
xenium_std <- FindNeighbors(xenium_std, dims = 1:50)
xenium_std <- FindClusters(xenium_std, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 36602
    ## Number of edges: 1016983
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9441
    ## Number of communities: 23
    ## Elapsed time: 3 seconds

``` r
n_clust <- length(levels(Idents(xenium_std)))
pal <- c(mypal, hcl.colors(max(0, n_clust - length(mypal))))
ImageDimPlot(xenium_std, size = 0.5, cols = pal) +
    coord_flip() + scale_x_reverse() + NoLegend() +
    theme(plot.background = element_rect(fill = 'white'),
          panel.background = element_rect(fill = 'white'))
```

<img src="banksy_files/figure-gfm/xenium_std_viz-1.png" alt="" style="display: block; margin: auto;" />

We compare the two embeddings with three checks:

``` r
U1 <- scale(Embeddings(xenium_std, 'pca'), scale = FALSE)
U2 <- scale(Embeddings(xenium, 'BANKSY'), scale = FALSE)
npcs <- ncol(U1)
```

**Subspace overlap.** Individual PCs may differ due to near-degenerate
singular values, but the subspaces they span should be equivalent. We
measure this via principal angles between the top-k subspaces from each
path. A minimum cosine near 1 means the two subspaces are effectively
identical:

``` r
Q1 <- qr.Q(qr(U1))
Q2 <- qr.Q(qr(U2))
subspace_cos <- sapply(seq_len(npcs), function(k) {
    min(svd(crossprod(Q1[, 1:k, drop = FALSE],
                      Q2[, 1:k, drop = FALSE]))$d)
})

plot(seq_len(npcs), subspace_cos, type = 'b', pch = 19, ylim = c(0.5, 1),
     xlab = 'Number of PCs (k)', ylab = 'Min cosine of principal angles',
     main = 'Subspace overlap: standard vs lazy')
abline(h = 0.99, lty = 2, col = 'grey50')
```

<img src="banksy_files/figure-gfm/xenium_subspace-1.png" alt="" style="display: block; margin: auto;" />

**kNN overlap.** Since clustering operates on the nearest-neighbor
graph, we compare the 30-NN graphs built from each embedding. The mean
Jaccard index measures what fraction of each cell’s neighbors are shared
between the two paths:

``` r
nn1 <- RANN::nn2(U1, k = 30)$nn.idx
nn2 <- RANN::nn2(U2, k = 30)$nn.idx
jaccard <- mean(sapply(seq_len(nrow(nn1)), function(i) {
    length(intersect(nn1[i,], nn2[i,])) / length(union(nn1[i,], nn2[i,]))
}))
cat('Mean kNN Jaccard overlap (k=30):', round(jaccard, 4), '\n')
```

    ## Mean kNN Jaccard overlap (k=30): 0.9999

**Clustering agreement.** The adjusted Rand index (ARI) measures how
closely the two sets of cluster labels agree, adjusted for chance:

``` r
ari <- mclust::adjustedRandIndex(
    as.integer(Idents(xenium_std)),
    as.integer(Idents(xenium))
)
cat('Adjusted Rand Index:', round(ari, 4), '\n')
```

    ## Adjusted Rand Index: 0.8503

</details>

## Getting help

For more information, visit <https://github.com/prabhakarlab/Banksy>.

<details>

<summary>

Vignette runtime
</summary>

    ## Time difference of 14.24326 mins

</details>

<details>

<summary>

Session info
</summary>

``` r
sessionInfo()
```

    ## R version 4.4.3 (2025-02-28)
    ## Platform: x86_64-conda-linux-gnu
    ## Running under: Rocky Linux 9.7 (Blue Onyx)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /gpfs/scrubbed/jxlee/tmp/codex/banksy-split-env-fresh/lib/libopenblasp-r0.3.33.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] harmony_2.0.2               Rcpp_1.1.1-1.1             
    ##  [3] ExperimentHub_2.14.0        AnnotationHub_3.14.0       
    ##  [5] BiocFileCache_2.14.0        dbplyr_2.5.2               
    ##  [7] spatialLIBD_1.18.0          SpatialExperiment_1.16.0   
    ##  [9] SingleCellExperiment_1.28.0 SummarizedExperiment_1.36.0
    ## [11] Biobase_2.66.0              GenomicRanges_1.58.0       
    ## [13] GenomeInfoDb_1.42.0         IRanges_2.40.0             
    ## [15] S4Vectors_0.44.0            BiocGenerics_0.52.0        
    ## [17] MatrixGenerics_1.18.0       matrixStats_1.5.0          
    ## [19] future_1.70.0               pals_1.10                  
    ## [21] gridExtra_2.3               ggplot2_4.0.3              
    ## [23] SeuratWrappers_0.4.0        ssHippo.SeuratData_3.1.4   
    ## [25] SeuratData_0.2.2.9002       Seurat_5.5.0               
    ## [27] SeuratObject_5.4.0          sp_2.2-1                   
    ## [29] Banksy_1.9.1               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] bitops_1.0-9             spatstat.sparse_3.2-0    httr_1.4.8              
    ##   [4] RColorBrewer_1.1-3       doParallel_1.0.17        tools_4.4.3             
    ##   [7] sctransform_0.4.3        R6_2.6.1                 DT_0.34.0               
    ##  [10] lazyeval_0.2.3           uwot_0.2.4               withr_3.0.2             
    ##  [13] progressr_0.19.0         cli_3.6.6                spatstat.explore_3.8-0  
    ##  [16] fastDummies_1.7.6        sass_0.4.10              labeling_0.4.3          
    ##  [19] arrow_24.0.0             S7_0.2.2                 spatstat.data_3.1-9     
    ##  [22] ggridges_0.5.7           pbapply_1.7-4            Rsamtools_2.22.0        
    ##  [25] dbscan_1.2.4             R.utils_2.13.0           aricode_1.1.0           
    ##  [28] scater_1.34.1            dichromat_2.0-0.1        sessioninfo_1.2.3       
    ##  [31] parallelly_1.47.0        attempt_0.3.1            maps_3.4.3              
    ##  [34] limma_3.62.1             RSQLite_3.53.1           BiocIO_1.16.0           
    ##  [37] generics_0.1.4           ica_1.0-3                spatstat.random_3.4-5   
    ##  [40] dplyr_1.2.1              Matrix_1.7-5             ggbeeswarm_0.7.3        
    ##  [43] abind_1.4-8              R.methodsS3_1.8.2        lifecycle_1.0.5         
    ##  [46] yaml_2.3.12              edgeR_4.4.0              SparseArray_1.6.0       
    ##  [49] Rtsne_0.17               paletteer_1.7.0          grid_4.4.3              
    ##  [52] blob_1.3.0               promises_1.5.0           crayon_1.5.3            
    ##  [55] miniUI_0.1.2             lattice_0.22-9           beachmat_2.22.0         
    ##  [58] cowplot_1.2.0            KEGGREST_1.46.0          mapproj_1.2.12          
    ##  [61] magick_2.9.1             pillar_1.11.1            knitr_1.51              
    ##  [64] rjson_0.2.23             future.apply_1.20.2      codetools_0.2-20        
    ##  [67] glue_1.8.1               spatstat.univar_3.2-0    data.table_1.17.8       
    ##  [70] remotes_2.5.0            vctrs_0.7.3              png_0.1-9               
    ##  [73] spam_2.11-3              gtable_0.3.6             assertthat_0.2.1        
    ##  [76] rematch2_2.1.2           cachem_1.1.0             xfun_0.57               
    ##  [79] S4Arrays_1.6.0           mime_0.13                survival_3.8-6          
    ##  [82] RcppHungarian_0.3        iterators_1.0.14         fields_17.3             
    ##  [85] statmod_1.5.2            fitdistrplus_1.2-6       ROCR_1.0-12             
    ##  [88] nlme_3.1-169             bit64_4.8.2              filelock_1.0.3          
    ##  [91] RcppAnnoy_0.0.23         bslib_0.11.0             irlba_2.3.7             
    ##  [94] vipor_0.4.7              KernSmooth_2.23-26       otel_0.2.0              
    ##  [97] colorspace_2.1-2         DBI_1.3.0                tidyselect_1.2.1        
    ## [100] bit_4.6.0                compiler_4.4.3           curl_7.1.0              
    ## [103] BiocNeighbors_2.0.0      hdf5r_1.3.12             DelayedArray_0.32.0     
    ## [106] plotly_4.12.0            rtracklayer_1.66.0       scales_1.4.0            
    ## [109] lmtest_0.9-40            rappdirs_0.3.4           stringr_1.6.0           
    ## [112] digest_0.6.39            goftest_1.2-3            spatstat.utils_3.2-3    
    ## [115] rmarkdown_2.31           benchmarkmeData_2.0.0    RhpcBLASctl_0.23-42     
    ## [118] XVector_0.46.0           htmltools_0.5.9          pkgconfig_2.0.3         
    ## [121] fastmap_1.2.0            rlang_1.2.0              htmlwidgets_1.6.4       
    ## [124] UCSC.utils_1.2.0         shiny_1.13.0             jquerylib_0.1.4         
    ## [127] farver_2.1.2             zoo_1.8-15               jsonlite_2.0.0          
    ## [130] BiocParallel_1.40.0      mclust_6.1.2             config_0.3.2            
    ## [133] R.oo_1.27.1              BiocSingular_1.22.0      RCurl_1.98-1.18         
    ## [136] magrittr_2.0.5           scuttle_1.16.0           GenomeInfoDbData_1.2.13 
    ## [139] dotCall64_1.2            patchwork_1.3.2          viridis_0.6.5           
    ## [142] reticulate_1.46.0        leidenAlg_1.1.6          stringi_1.8.7           
    ## [145] zlibbioc_1.52.0          MASS_7.3-65              plyr_1.8.9              
    ## [148] parallel_4.4.3           listenv_0.10.1           ggrepel_0.9.8           
    ## [151] deldir_2.0-4             Biostrings_2.74.0        sccore_1.0.7            
    ## [154] splines_4.4.3            tensor_1.5.1             locfit_1.5-9.12         
    ## [157] igraph_2.3.1             spatstat.geom_3.8-1      RcppHNSW_0.6.0          
    ## [160] ScaledMatrix_1.14.0      reshape2_1.4.5           BiocVersion_3.20.0      
    ## [163] XML_3.99-0.23            evaluate_1.0.5           golem_0.5.1             
    ## [166] BiocManager_1.30.27      foreach_1.5.2            httpuv_1.6.17           
    ## [169] RANN_2.6.2               tidyr_1.3.2              purrr_1.2.2             
    ## [172] polyclip_1.10-7          benchmarkme_1.0.8        scattermore_1.2         
    ## [175] rsvd_1.0.5               xtable_1.8-8             restfulr_0.0.16         
    ## [178] RSpectra_0.16-2          later_1.4.8              viridisLite_0.4.3       
    ## [181] tibble_3.3.1             beeswarm_0.4.0           GenomicAlignments_1.42.0
    ## [184] memoise_2.0.1            AnnotationDbi_1.68.0     cluster_2.1.8.2         
    ## [187] shinyWidgets_0.9.1       globals_0.19.1

</details>
