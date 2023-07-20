Running BANKSY with Seurat
================
Compiled: July 20, 2023

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#overview" id="toc-overview">Overview</a>
-   <a href="#running-banksy-within-seurats-spatial-framework"
    id="toc-running-banksy-within-seurats-spatial-framework">Running BANKSY
    within Seurat’s spatial framework</a>
-   <a href="#running-banksy-with-locations-provided-explicitly"
    id="toc-running-banksy-with-locations-provided-explicitly">Running
    BANKSY with locations provided explicitly</a>
-   <a href="#getting-help" id="toc-getting-help">Getting help</a>

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

**Caveat**: `ScaleData` should not be run after a call to `RunBanksy`;
`RunBanksy` populates the `scale.data` slot with the scaled BANKSY
matrix. Calling `ScaleData` after `RunBanksy` performs gene-wise
z-scaling, negating the effect of `lambda`.

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
ss.hippo <- RunBanksy(ss.hippo, lambda = 0.2, verbose=TRUE, 
                      assay = 'Spatial', slot = 'data', features = 'variable',
                      k_geom = 15)
ss.hippo
```

    ## An object of class Seurat 
    ## 27264 features across 10000 samples within 2 assays 
    ## Active assay: BANKSY (4000 features, 0 variable features)
    ##  1 other assay present: Spatial
    ##  1 image present: image

Note that the `RunBanksy` function sets the default assay to `BANKSY` (
determined by the `assay_name` argument) and fills the `scale.data`
slot.

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
    ## Number of edges: 366094
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9035
    ## Number of communities: 13
    ## Elapsed time: 2 seconds

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
    ## SNAP25 3.381910e-46  -1.207905 0.657 0.824 7.867674e-42
    ## CHGB   1.791703e-44  -1.873680 0.437 0.698 4.168217e-40
    ## STMN2  6.651615e-25  -1.263710 0.334 0.576 1.547432e-20
    ## SYN2   6.404044e-23  -1.380049 0.334 0.566 1.489837e-18
    ## ATP2B1 1.298411e-21   1.150366 0.636 0.479 3.020623e-17
    ## CPLX2  3.698909e-21  -1.048808 0.287 0.523 8.605141e-17
    ## PCP4   4.083264e-18  -1.058147 0.379 0.578 9.499305e-14
    ## PRKCB  5.548792e-18   1.224930 0.547 0.341 1.290871e-13
    ## DDN    4.230865e-14   1.290755 0.589 0.396 9.842685e-10

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

Scale the BANKSY matrix and run PCA:

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

    ## Time difference of 2.153743 mins

</details>
<details>
<summary>
Session info
</summary>

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19043)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Singapore.utf8  LC_CTYPE=English_Singapore.utf8   
    ## [3] LC_MONETARY=English_Singapore.utf8 LC_NUMERIC=C                      
    ## [5] LC_TIME=English_Singapore.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pals_1.7                 gridExtra_2.3            ggplot2_3.4.0           
    ##  [4] SeuratWrappers_0.3.0     ssHippo.SeuratData_3.1.4 pbmc3k.SeuratData_3.1.4 
    ##  [7] SeuratData_0.2.2         SeuratObject_4.1.3       Seurat_4.3.0            
    ## [10] Banksy_0.1.5            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2                  spatstat.explore_3.0-5     
    ##   [3] reticulate_1.27             R.utils_2.12.2             
    ##   [5] tidyselect_1.2.0            htmlwidgets_1.6.0          
    ##   [7] grid_4.2.1                  Rtsne_0.16                 
    ##   [9] munsell_0.5.0               codetools_0.2-18           
    ##  [11] ica_1.0-3                   future_1.30.0              
    ##  [13] miniUI_0.1.1.1              withr_2.5.0                
    ##  [15] spatstat.random_3.0-1       colorspace_2.0-3           
    ##  [17] progressr_0.12.0            Biobase_2.58.0             
    ##  [19] highr_0.10                  knitr_1.41                 
    ##  [21] ggalluvial_0.12.3           rstudioapi_0.14            
    ##  [23] stats4_4.2.1                ROCR_1.0-11                
    ##  [25] tensor_1.5                  listenv_0.9.0              
    ##  [27] MatrixGenerics_1.10.0       labeling_0.4.2             
    ##  [29] GenomeInfoDbData_1.2.9      polyclip_1.10-4            
    ##  [31] farver_2.1.1                parallelly_1.33.0          
    ##  [33] vctrs_0.5.1                 generics_0.1.3             
    ##  [35] xfun_0.36                   R6_2.5.1                   
    ##  [37] doParallel_1.0.17           GenomeInfoDb_1.34.4        
    ##  [39] clue_0.3-62                 rsvd_1.0.5                 
    ##  [41] bitops_1.0-7                spatstat.utils_3.0-1       
    ##  [43] DelayedArray_0.24.0         assertthat_0.2.1           
    ##  [45] promises_1.2.0.1            scales_1.2.1               
    ##  [47] gtable_0.3.1                globals_0.16.2             
    ##  [49] goftest_1.2-3               rlang_1.1.1                
    ##  [51] GlobalOptions_0.1.2         splines_4.2.1              
    ##  [53] lazyeval_0.2.2              dichromat_2.0-0.1          
    ##  [55] spatstat.geom_3.0-3         BiocManager_1.30.19        
    ##  [57] yaml_2.3.6                  reshape2_1.4.4             
    ##  [59] abind_1.4-5                 httpuv_1.6.6               
    ##  [61] tools_4.2.1                 sccore_1.0.2               
    ##  [63] ellipsis_0.3.2              RColorBrewer_1.1-3         
    ##  [65] BiocGenerics_0.44.0         ggridges_0.5.4             
    ##  [67] Rcpp_1.0.9                  plyr_1.8.8                 
    ##  [69] progress_1.2.2              zlibbioc_1.44.0            
    ##  [71] purrr_0.3.5                 RCurl_1.98-1.9             
    ##  [73] prettyunits_1.1.1           dbscan_1.1-11              
    ##  [75] deldir_1.0-6                pbapply_1.7-0              
    ##  [77] GetoptLong_1.0.5            cowplot_1.1.1              
    ##  [79] S4Vectors_0.36.0            zoo_1.8-11                 
    ##  [81] SummarizedExperiment_1.28.0 ggrepel_0.9.2              
    ##  [83] cluster_2.1.4               magrittr_2.0.3             
    ##  [85] data.table_1.14.6           scattermore_0.8            
    ##  [87] circlize_0.4.15             lmtest_0.9-40              
    ##  [89] RANN_2.6.1                  fitdistrplus_1.1-8         
    ##  [91] matrixStats_0.62.0          hms_1.1.2                  
    ##  [93] patchwork_1.1.2             mime_0.12                  
    ##  [95] evaluate_0.19               xtable_1.8-4               
    ##  [97] mclust_6.0.0                IRanges_2.32.0             
    ##  [99] shape_1.4.6                 compiler_4.2.1             
    ## [101] tibble_3.1.8                maps_3.4.1                 
    ## [103] KernSmooth_2.23-20          crayon_1.5.2               
    ## [105] R.oo_1.25.0                 htmltools_0.5.4            
    ## [107] later_1.3.0                 tidyr_1.2.1                
    ## [109] DBI_1.1.3                   ComplexHeatmap_2.14.0      
    ## [111] MASS_7.3-58.1               rappdirs_0.3.3             
    ## [113] Matrix_1.5-3                cli_3.4.1                  
    ## [115] R.methodsS3_1.8.2           parallel_4.2.1             
    ## [117] RcppHungarian_0.2           igraph_1.3.5               
    ## [119] GenomicRanges_1.50.1        pkgconfig_2.0.3            
    ## [121] sp_1.5-1                    plotly_4.10.1              
    ## [123] spatstat.sparse_3.0-0       foreach_1.5.2              
    ## [125] XVector_0.38.0              leidenAlg_1.1.0            
    ## [127] stringr_1.5.0               digest_0.6.30              
    ## [129] sctransform_0.3.5           RcppAnnoy_0.0.20           
    ## [131] spatstat.data_3.0-0         rmarkdown_2.19             
    ## [133] leiden_0.4.3                uwot_0.1.14                
    ## [135] shiny_1.7.4                 rjson_0.2.21               
    ## [137] lifecycle_1.0.3             nlme_3.1-160               
    ## [139] jsonlite_1.8.3              mapproj_1.2.9              
    ## [141] viridisLite_0.4.1           limma_3.54.0               
    ## [143] fansi_1.0.3                 pillar_1.8.1               
    ## [145] lattice_0.20-45             fastmap_1.1.0              
    ## [147] httr_1.4.4                  survival_3.4-0             
    ## [149] glue_1.6.2                  remotes_2.4.2              
    ## [151] png_0.1-7                   iterators_1.0.14           
    ## [153] stringi_1.7.8               dplyr_1.0.10               
    ## [155] irlba_2.3.5.1               future.apply_1.10.0

</details>
