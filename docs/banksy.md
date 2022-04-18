Running BANKSY with Seurat
================
Compiled: April 18, 2022

-   [Introduction](#introduction)
-   [Overview](#overview)
-   [Running BANKSY within Seurat’s spatial
    framework](#running-banksy-within-seurats-spatial-framework)
-   [Running BANKSY with locations provided
    explicitly](#running-banksy-with-locations-provided-explicitly)
-   [Getting help](#getting-help)

## Introduction

In this vignette, we describe how to run BANKSY with Seurat objects. If
you use BANKSY in your research, please cite

> *BANKSY: A Spatial Omics Algorithm that Unifies Cell Type Clustering
> and Tissue Domain Segmentation*
>
> Vipul Singhal, Nigel Chou, Joseph Lee, Jinyue Liu, Wan Kee Chock, Li
> Lin, Yun-Ching Chang, Erica Teo, Hwee Kuan Lee, Kok Hao Chen, Shyam
> Prabhakar
>
> bioRxiv, 2022
>
> doi:
> [10.1101/2022.04.14.488259](https://www.biorxiv.org/content/10.1101/2022.04.14.488259)
>
> Website: <https://prabhakarlab.github.io/Banksy>

BANKSY is a method that incorporates neighborhood information for
clustering spatial omics data. By doing so, BANKSY is able to

-   improve cell-type assignment in noisy data
-   distinguish subtly different cell-types stratified by
    microenvironment
-   identify spatial domains sharing the same microenvironment

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

Prerequisites to install:

-   [Seurat](https://satijalab.org/seurat/install)
-   [SeuratData](https://github.com/satijalab/seurat-data)
-   [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
-   [Banksy](https://github.com/prabhakarlab/Banksy/)

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
```

To run BANKSY, we specify the following:

-   `lambda`: a numeric value in \[0,1\]. With low values of lambda,
    BANKSY operates in cell-typing mode, while high values of lambda
    find spatial domains.
-   `assay` and `slot`: determines where to pull the expression data
    from
-   `features`: specifies features for downstream analysis. This can be
    `'all'`, `'variable'` or a subset of features.  
-   `k_geom`: the number of neighbors that defines a cell’s neighborhood

Call `?RunBanksy` for more details on function parameters.

``` r
# Run BANKSY
ss.hippo <- RunBanksy(ss.hippo, lambda = 0.15, verbose=TRUE, 
                      assay = 'Spatial', slot = 'data', features = 'variable',
                      k_geom = 10)
ss.hippo
```

    ## An object of class Seurat 
    ## 27264 features across 10000 samples within 2 assays 
    ## Active assay: BANKSY (4000 features, 0 variable features)
    ##  1 other assay present: Spatial
    ##  1 image present: image

Note that the `RunBanksy` function sets the default assay to `BANKSY` (
determined by the `assay_name` argument).

The rest of the pipeline is similar to the ‘default’ Seurat pipline. We
scale the data and run dimensionality reduction with PCA and UMAP:

``` r
# Scale
ss.hippo <- ScaleData(ss.hippo)
# Run PCA and UMAP
ss.hippo <- RunPCA(ss.hippo, assay = 'BANKSY', features = rownames(ss.hippo), npcs = 50)
ss.hippo <- RunUMAP(ss.hippo, dims = 1:30)
```

Next, find BANKSY clusters:

``` r
# Clustering
ss.hippo <- FindNeighbors(ss.hippo, dims = 1:30)
ss.hippo <- FindClusters(ss.hippo, resolution = 0.6)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10000
    ## Number of edges: 330767
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8985
    ## Number of communities: 15
    ## Elapsed time: 2 seconds

Visualize the UMAP and Spatial plot:

``` r
# Viz
grid.arrange(
    DimPlot(ss.hippo, cols = mypal, pt.size = 0.25),
    SpatialDimPlot(ss.hippo, stroke = NA, label = TRUE, label.size = 3, 
                   repel = TRUE, alpha = 0.5, cols = mypal, pt.size.factor = 3),
    ncol = 2
)
```

<img src="banksy_files/figure-gfm/ss_viz-1.png" style="display: block; margin: auto;" />

Find markers based on the BANKSY clusters and visualize them. Here, we
find differentially expressed genes between the CA1 and CA3 regions.

``` r
# Find markers
DefaultAssay(ss.hippo) <- 'Spatial'
markers <- FindMarkers(ss.hippo, ident.1 = 4, ident.2 = 5, only.pos = F, 
                       logfc.threshold = 1, min.pct = 0.5)
markers <- markers[markers$p_val_adj < 0.01,]
markers
```

    ##               p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## TMSB4X 9.721010e-57   1.209878 0.870 0.688 2.261496e-52
    ## SNAP25 7.547825e-28  -1.003154 0.656 0.748 1.755926e-23
    ## ATP2B1 3.799825e-26   1.166903 0.631 0.423 8.839913e-22
    ## CHGB   7.225647e-26  -1.701535 0.438 0.599 1.680975e-21
    ## SYN2   1.529839e-20  -1.370905 0.339 0.540 3.559018e-16
    ## DDN    1.688286e-20   1.344971 0.595 0.349 3.927629e-16
    ## PRKCB  1.067433e-18   1.196375 0.545 0.336 2.483276e-14
    ## WIPF3  5.379919e-18   1.142872 0.516 0.297 1.251584e-13
    ## STMN2  1.720493e-15  -1.112610 0.341 0.503 4.002554e-11

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
vf.hippo <- RunBanksy(vf.hippo, lambda = 0.3, dimx = 'sdimx', dimy = 'sdimy', 
                      assay = 'RNA', slot = 'data', features = 'all', k_geom = 10)
```

Scale the BANKSY matrix and run PCA:

``` r
# Scale
vf.hippo <- ScaleData(vf.hippo)
# PCA
vf.hippo <- RunPCA(vf.hippo, assay = 'BANKSY', features = rownames(vf.hippo), npcs = 50)
```

Find BANKSY clusters:

``` r
# Cluster
vf.hippo <- FindNeighbors(vf.hippo, dims = 1:30)
vf.hippo <- FindClusters(vf.hippo, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10205
    ## Number of edges: 410260
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8957
    ## Number of communities: 16
    ## Elapsed time: 2 seconds

Visualise BANKSY clusters in spatial dimensions:

``` r
# Viz
FeatureScatter(vf.hippo, 'sdimx', 'sdimy', cols = mypal, pt.size = 0.75)
```

<img src="banksy_files/figure-gfm/hippo_viz-1.png" style="display: block; margin: auto;" />

Find markers and visualise them. Here, we do so for a cluster defined by
a thin layer of cells expressing Gfap. We also write a simple function
`genePlot` that plots marker genes in spatial dimensions.

``` r
# Find markers
DefaultAssay(vf.hippo) <- 'RNA'
markers <- FindMarkers(vf.hippo, ident.1 = 6, only.pos = TRUE)

genePlot <- function(object, dimx, dimy, gene, 
                     slot = 'scale.data', q.low = 0.01, q.high = 0.99,
                     col.low='blue', col.high='red') {
    val <- GetAssayData(object, slot)[gene,]
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

## Getting help

For more information, visit <https://github.com/prabhakarlab/Banksy>.

<details>
<summary>
Vignette runtime
</summary>

    ## Time difference of 2.68403 mins

</details>
<details>
<summary>
Session info
</summary>

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19043)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Singapore.1252  LC_CTYPE=English_Singapore.1252    LC_MONETARY=English_Singapore.1252
    ## [4] LC_NUMERIC=C                       LC_TIME=English_Singapore.1252    
    ## system code page: 65001
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pals_1.7                  gridExtra_2.3             ggplot2_3.3.5             SeuratWrappers_0.3.0     
    ##  [5] stxBrain.SeuratData_0.1.1 ssHippo.SeuratData_3.1.4  SeuratData_0.2.1          sp_1.4-6                 
    ##  [9] SeuratObject_4.0.999.9011 Seurat_4.1.0.9004         Banksy_0.1.3             
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2                  reticulate_1.24             R.utils_2.11.0             
    ##   [4] tidyselect_1.1.2            htmlwidgets_1.5.4           grid_4.1.2                 
    ##   [7] Rtsne_0.15                  munsell_0.5.0               codetools_0.2-18           
    ##  [10] ica_1.0-2                   future_1.24.0               miniUI_0.1.1.1             
    ##  [13] withr_2.5.0                 spatstat.random_2.1-0       colorspace_2.0-2           
    ##  [16] progressr_0.10.0            Biobase_2.54.0              highr_0.9                  
    ##  [19] knitr_1.37                  ggalluvial_0.12.3           rstudioapi_0.13            
    ##  [22] stats4_4.1.2                ROCR_1.0-11                 tensor_1.5                 
    ##  [25] listenv_0.8.0               labeling_0.4.2              MatrixGenerics_1.6.0       
    ##  [28] GenomeInfoDbData_1.2.7      polyclip_1.10-0             farver_2.1.0               
    ##  [31] parallelly_1.30.0           Matrix.utils_0.9.8          vctrs_0.3.8                
    ##  [34] generics_0.1.2              xfun_0.29                   R6_2.5.1                   
    ##  [37] doParallel_1.0.17           GenomeInfoDb_1.30.1         clue_0.3-60                
    ##  [40] rsvd_1.0.5                  bitops_1.0-7                spatstat.utils_2.3-0       
    ##  [43] DelayedArray_0.20.0         assertthat_0.2.1            promises_1.2.0.1           
    ##  [46] scales_1.1.1                rgeos_0.5-9                 gtable_0.3.0               
    ##  [49] globals_0.14.0              goftest_1.2-3               rlang_1.0.2                
    ##  [52] GlobalOptions_0.1.2         splines_4.1.2               lazyeval_0.2.2             
    ##  [55] dichromat_2.0-0             spatstat.geom_2.3-2         BiocManager_1.30.16        
    ##  [58] yaml_2.2.1                  reshape2_1.4.4              abind_1.4-5                
    ##  [61] httpuv_1.6.5                tools_4.1.2                 sccore_1.0.1               
    ##  [64] ellipsis_0.3.2              spatstat.core_2.4-0         RColorBrewer_1.1-2         
    ##  [67] BiocGenerics_0.40.0         ggridges_0.5.3              Rcpp_1.0.7                 
    ##  [70] plyr_1.8.6                  zlibbioc_1.40.0             purrr_0.3.4                
    ##  [73] RCurl_1.98-1.6              rpart_4.1-15                dbscan_1.1-10              
    ##  [76] deldir_1.0-6                pbapply_1.5-0               GetoptLong_1.0.5           
    ##  [79] cowplot_1.1.1               S4Vectors_0.32.3            zoo_1.8-9                  
    ##  [82] SummarizedExperiment_1.24.0 grr_0.9.5                   ggrepel_0.9.1              
    ##  [85] cluster_2.1.2               magrittr_2.0.1              RSpectra_0.16-0            
    ##  [88] data.table_1.14.2           scattermore_0.8             circlize_0.4.14            
    ##  [91] lmtest_0.9-39               RANN_2.6.1                  fitdistrplus_1.1-8         
    ##  [94] matrixStats_0.61.0          patchwork_1.1.1             mime_0.12                  
    ##  [97] evaluate_0.15               xtable_1.8-4                mclust_5.4.9               
    ## [100] IRanges_2.28.0              shape_1.4.6                 compiler_4.1.2             
    ## [103] tibble_3.1.6                maps_3.4.0                  KernSmooth_2.23-20         
    ## [106] crayon_1.5.0                R.oo_1.24.0                 htmltools_0.5.2            
    ## [109] mgcv_1.8-39                 later_1.3.0                 tidyr_1.2.0                
    ## [112] DBI_1.1.2                   ComplexHeatmap_2.10.0       MASS_7.3-54                
    ## [115] rappdirs_0.3.3              Matrix_1.3-4                cli_3.1.0                  
    ## [118] R.methodsS3_1.8.1           parallel_4.1.2              RcppHungarian_0.2          
    ## [121] igraph_1.2.11               GenomicRanges_1.46.1        pkgconfig_2.0.3            
    ## [124] plotly_4.10.0               spatstat.sparse_2.1-0       foreach_1.5.2              
    ## [127] XVector_0.34.0              leidenAlg_1.0.2             stringr_1.4.0              
    ## [130] digest_0.6.29               sctransform_0.3.3           RcppAnnoy_0.0.19           
    ## [133] spatstat.data_2.1-2         rmarkdown_2.13              leiden_0.3.9               
    ## [136] uwot_0.1.11                 shiny_1.7.1                 rjson_0.2.21               
    ## [139] lifecycle_1.0.1             nlme_3.1-153                jsonlite_1.8.0             
    ## [142] mapproj_1.2.8               limma_3.50.1                viridisLite_0.4.0          
    ## [145] fansi_0.5.0                 pillar_1.7.0                lattice_0.20-45            
    ## [148] fastmap_1.1.0               httr_1.4.2                  survival_3.2-13            
    ## [151] glue_1.6.0                  remotes_2.4.2               png_0.1-7                  
    ## [154] iterators_1.0.14            stringi_1.7.6               dplyr_1.0.7                
    ## [157] irlba_2.3.5                 future.apply_1.8.1

</details>
