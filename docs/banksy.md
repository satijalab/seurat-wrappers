Running BANKSY with Seurat
================
Compiled: April 02, 2024

- [Introduction](#introduction)
- [Overview](#overview)
- [Running BANKSY within Seurat’s spatial
  framework](#running-banksy-within-seurats-spatial-framework)
- [Running BANKSY with locations provided
  explicitly](#running-banksy-with-locations-provided-explicitly)
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
    ## Number of edges: 366882
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9034
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
    ## SNAP25 3.892818e-45  -1.232135 0.660 0.820 9.056251e-41
    ## CHGB   2.284338e-43  -1.981091 0.439 0.694 5.314283e-39
    ## STMN2  1.579733e-24  -1.432890 0.334 0.572 3.675092e-20
    ## ATP2B1 8.175474e-23   1.253406 0.639 0.472 1.901942e-18
    ## SYN2   1.532794e-22  -1.557895 0.334 0.562 3.565892e-18
    ## CPLX2  1.130974e-20  -1.202986 0.289 0.520 2.631099e-16
    ## PRKCB  2.980819e-18   1.386395 0.548 0.340 6.934577e-14
    ## PCP4   4.140456e-18  -1.194483 0.381 0.578 9.632356e-14
    ## TUBB2A 5.408558e-16  -1.044764 0.448 0.622 1.258247e-11
    ## DDN    3.016969e-14   1.398636 0.589 0.394 7.018677e-10
    ## SNCA   1.301494e-11  -1.012388 0.397 0.542 3.027796e-07

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
    ## Number of edges: 447627
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9108
    ## Number of communities: 15
    ## Elapsed time: 2 seconds

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

## Getting help

For more information, visit <https://github.com/prabhakarlab/Banksy>.

<details>
<summary>
Vignette runtime
</summary>

    ## Time difference of 1.316144 mins

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
    ## [10] Banksy_0.1.6            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppHungarian_0.3           RcppAnnoy_0.0.22           
    ##   [3] splines_4.3.2               later_1.3.2                
    ##   [5] bitops_1.0-7                R.oo_1.26.0                
    ##   [7] tibble_3.2.1                polyclip_1.10-6            
    ##   [9] fastDummies_1.7.3           lifecycle_1.0.4            
    ##  [11] doParallel_1.0.17           globals_0.16.2             
    ##  [13] lattice_0.22-5              MASS_7.3-60.0.1            
    ##  [15] magrittr_2.0.3              limma_3.58.1               
    ##  [17] plotly_4.10.4               rmarkdown_2.25             
    ##  [19] remotes_2.4.2.1             yaml_2.3.8                 
    ##  [21] httpuv_1.6.14               sctransform_0.4.1          
    ##  [23] spam_2.10-0                 spatstat.sparse_3.0-3      
    ##  [25] reticulate_1.35.0           cowplot_1.1.3              
    ##  [27] mapproj_1.2.11              pbapply_1.7-2              
    ##  [29] RColorBrewer_1.1-3          maps_3.4.2                 
    ##  [31] abind_1.4-5                 zlibbioc_1.48.0            
    ##  [33] Rtsne_0.17                  GenomicRanges_1.54.1       
    ##  [35] R.utils_2.12.3              purrr_1.0.2                
    ##  [37] BiocGenerics_0.48.1         RCurl_1.98-1.14            
    ##  [39] rappdirs_0.3.3              circlize_0.4.15            
    ##  [41] GenomeInfoDbData_1.2.11     IRanges_2.36.0             
    ##  [43] S4Vectors_0.40.2            ggrepel_0.9.5              
    ##  [45] irlba_2.3.5.1               listenv_0.9.1              
    ##  [47] spatstat.utils_3.0-4        goftest_1.2-3              
    ##  [49] RSpectra_0.16-1             spatstat.random_3.2-2      
    ##  [51] fitdistrplus_1.1-11         parallelly_1.37.0          
    ##  [53] leiden_0.4.3.1              codetools_0.2-19           
    ##  [55] DelayedArray_0.28.0         tidyselect_1.2.0           
    ##  [57] shape_1.4.6                 farver_2.1.1               
    ##  [59] matrixStats_1.2.0           stats4_4.3.2               
    ##  [61] spatstat.explore_3.2-6      jsonlite_1.8.8             
    ##  [63] GetoptLong_1.0.5            ellipsis_0.3.2             
    ##  [65] progressr_0.14.0            ggridges_0.5.6             
    ##  [67] ggalluvial_0.12.5           survival_3.5-7             
    ##  [69] iterators_1.0.14            foreach_1.5.2              
    ##  [71] dbscan_1.1-12               tools_4.3.2                
    ##  [73] progress_1.2.3              ica_1.0-3                  
    ##  [75] Rcpp_1.0.12                 glue_1.7.0                 
    ##  [77] SparseArray_1.2.4           xfun_0.42                  
    ##  [79] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.6        
    ##  [81] dplyr_1.1.4                 withr_3.0.0                
    ##  [83] BiocManager_1.30.22         fastmap_1.1.1              
    ##  [85] fansi_1.0.6                 rsvd_1.0.5                 
    ##  [87] digest_0.6.34               R6_2.5.1                   
    ##  [89] mime_0.12                   colorspace_2.1-0           
    ##  [91] scattermore_1.2             sccore_1.0.4               
    ##  [93] tensor_1.5                  dichromat_2.0-0.1          
    ##  [95] spatstat.data_3.0-4         R.methodsS3_1.8.2          
    ##  [97] utf8_1.2.4                  tidyr_1.3.1                
    ##  [99] generics_0.1.3              data.table_1.15.0          
    ## [101] prettyunits_1.2.0           httr_1.4.7                 
    ## [103] htmlwidgets_1.6.4           S4Arrays_1.2.0             
    ## [105] uwot_0.1.16                 pkgconfig_2.0.3            
    ## [107] gtable_0.3.4                ComplexHeatmap_2.18.0      
    ## [109] lmtest_0.9-40               XVector_0.42.0             
    ## [111] htmltools_0.5.7             dotCall64_1.1-1            
    ## [113] clue_0.3-65                 scales_1.3.0               
    ## [115] Biobase_2.62.0              png_0.1-8                  
    ## [117] knitr_1.45                  rstudioapi_0.15.0          
    ## [119] reshape2_1.4.4              rjson_0.2.21               
    ## [121] nlme_3.1-164                zoo_1.8-12                 
    ## [123] GlobalOptions_0.1.2         stringr_1.5.1              
    ## [125] KernSmooth_2.23-22          parallel_4.3.2             
    ## [127] miniUI_0.1.1.1              pillar_1.9.0               
    ## [129] grid_4.3.2                  vctrs_0.6.5                
    ## [131] RANN_2.6.1                  promises_1.2.1             
    ## [133] xtable_1.8-4                cluster_2.1.6              
    ## [135] evaluate_0.23               cli_3.6.2                  
    ## [137] compiler_4.3.2              rlang_1.1.3                
    ## [139] crayon_1.5.2                future.apply_1.11.1        
    ## [141] labeling_0.4.3              mclust_6.0.1               
    ## [143] plyr_1.8.9                  stringi_1.8.3              
    ## [145] viridisLite_0.4.2           deldir_2.0-2               
    ## [147] munsell_0.5.0               lazyeval_0.2.2             
    ## [149] spatstat.geom_3.2-8         Matrix_1.6-5               
    ## [151] RcppHNSW_0.6.0              hms_1.1.3                  
    ## [153] patchwork_1.2.0             future_1.33.1              
    ## [155] statmod_1.5.0               shiny_1.8.0                
    ## [157] highr_0.10                  SummarizedExperiment_1.32.0
    ## [159] ROCR_1.0-11                 leidenAlg_1.1.2            
    ## [161] igraph_2.0.1.1

</details>
