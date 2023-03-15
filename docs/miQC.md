Running miQC on Seurat objects
================
Compiled: July 19, 2021

This vigettte demonstrates the use of the miQC package in Seurat.
Vignette is based off of the [miQC
vignette](https://github.com/greenelab/miQC). If you use miQC in your
work, please cite:

> *miQC: An adaptive probabilistic framework for quality control of
> single-cell RNA-sequencing data*
> 
> Ariel A. Hippen, Matias M. Falco, Lukas M. Weber, Erdogan Pekcan
> Erkan, Kaiyang Zhang, Jennifer Anne Doherty, Anna Vähärautio, Casey S.
> Greene, Stephanie C. Hicks
> 
> bioRxiv, 2021
> 
> doi:
> [10.1101/2021.03.03.433798](https://www.biorxiv.org/content/10.1101/2021.03.03.433798v1)
> 
> GitHub: <https://github.com/greenelab/miQC>

Prerequisites to install:

  - [Seurat](https://satijalab.org/seurat/install)
  - [SeuratData](https://github.com/satijalab/seurat-data)
  - [flexmix](https://cran.r-project.org/web/packages/flexmix/index.html)
    which is wrapped by the [miQC](https://github.com/greenelab/miQC)
    package.
      - *At this point, the miQC algorithm has been adapted for use in
        Seurat through installation of flexmix only*.

<!-- end list -->

``` r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(flexmix)
```

## Introduction

This vignette provides a basic example of how to run miQC, which allows
users to perform cell-wise filtering of single-cell RNA-seq data for
quality control. Single-cell RNA-seq data is very sensitive to tissue
quality and choice of experimental workflow; it’s critical to ensure
compromised cells and failed cell libraries are removed. A high
proportion of reads mapping to mitochondrial DNA is one sign of a
damaged cell, so most analyses will remove cells with mtRNA over a
certain threshold, but those thresholds can be arbitrary and/or
detrimentally stringent, especially for archived tumor tissues. miQC
jointly models both the proportion of reads mapping to mtDNA genes and
the number of detected genes with mixture models in a probabilistic
framework to identify the low-quality cells in a given dataset.

## Example data

To demonstrate how to run miQC on a single-cell RNA-seq dataset, we’ll
use the `pbmc3k`dataset from the SeuratData package.

``` r
InstallData("pbmc3k")
data("pbmc3k")
pbmc3k
```

    ## An object of class Seurat 
    ## 13714 features across 2700 samples within 1 assay 
    ## Active assay: RNA (13714 features, 0 variable features)

## Seurat preprocessing

*miQC* requires two QC metrics for each single cell dataset: (1) the
number of unique genes detected per cell and (2) the percent
mitochondrial reads. The number of unique genes detected per cell are
typically calculated and stored automatically as metadata
(*nFeature\_RNA*) upon creation of a Seurat object with
`CreateSeuratObject`.

In order to calculate the percent mitochondrial reads in a cell we can
use `PercentageFeatureSet`. Human mitochondrial genes start with *MT-*
(and *mt-* for murine genes). For other IDs, we recommend using a
*biomaRt* query to map to chromosomal location and identify all
mitochondrial genes. We add this as metadata here to the Seurat object
as `"percent.mt"`.

``` r
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(object = pbmc3k, pattern = "^MT-")
```

## miQC

We can visually inspect the `"percent.mt"` and `"nFeature_RNA"` values
in the `pbmc3k` dataset.

``` r
FeatureScatter(pbmc3k, feature1 = "nFeature_RNA", feature2 = "percent.mt")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/miQC_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

We can see that most cells have a fairly low proportion of mitochondrial
reads, given that the graph is much denser at the bottom. We likely have
many cells that are intact and biologically meaningful. There are also a
few cells that have almost half of their reads mapping to mitochondrial
genes, which are likely broken or otherwise compromised and we will want
to exclude from our downstream analysis. However, it’s not clear what
boundaries to draw to separate the two groups of cells. With that in
mind, we’ll generate a linear mixture model using the `RunMiQC`
function. The linear mixture model will be stored in the `misc` slot of
the Seurat object as `"flexmix_model"`.

``` r
pbmc3k <- RunMiQC(pbmc3k, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75, 
    model.slot = "flexmix_model")
```

This function is a wrapper for *flexmix*, which fits a mixture model on
our data and returns the parameters of the two lines that best fit the
data, as well as the posterior probability of each cell being derived
from each distribution.

We can look at the parameters and posterior values directly with the
functions

``` r
flexmix::parameters(Misc(pbmc3k, "flexmix_model"))
```

    ##                         Comp.1       Comp.2
    ## coef.(Intercept)  2.005173e+00  7.144404951
    ## coef.nFeature_RNA 3.205993e-05 -0.004140063
    ## sigma             7.410023e-01  2.122227100

``` r
head(flexmix::posterior(Misc(pbmc3k, "flexmix_model")))
```

    ##           [,1]       [,2]
    ## [1,] 0.9288831 0.07111692
    ## [2,] 0.7604008 0.23959917
    ## [3,] 0.9196259 0.08037411
    ## [4,] 0.9711273 0.02887266
    ## [5,] 0.9873905 0.01260947
    ## [6,] 0.9782491 0.02175086

Or we can visualize the model results using the *PlotMiQC* function,
where `"miQC.probability"` represents the posterior probability of the
cell belonging to the compromised condition:

``` r
PlotMiQC(pbmc3k, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/miQC_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

As expected, the cells at the very top of the graph are almost certainly
compromised, most likely to have been derived from the distribution with
fewer unique genes and higher baseline mitochondrial expression.

We can use these posterior probabilities to choose which cells to keep,
and visualize the consequences of this filtering with the *PlotMiQC*
function. Recall when running `"RunMiQC"` we set the
`"posterior.cutoff"` to be 0.75.

``` r
PlotMiQC(pbmc3k, color.by = "miQC.keep")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/miQC_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

To actually perform the filtering and remove the indicated cells from
our Seurat object, we can subset the Seurat object parameter as such:

``` r
pbmc3k_filtered <- subset(pbmc3k, miQC.keep == "keep")
pbmc3k_filtered
```

    ## An object of class Seurat 
    ## 13714 features across 2593 samples within 1 assay 
    ## Active assay: RNA (13714 features, 0 variable features)

## Extras

In most cases, a linear mixture model will be satisfactory as well as
simplest, but *RunMiQC* also supports some non-linear mixture models:
currently polynomials and b-splines. A user should only need to change
the *model.type* parameter when making the model, and all visualization
and filtering functions will work the same as with a linear model.

``` r
pbmc3k <- RunMiQC(pbmc3k, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75, 
    model.slot = "flexmix_model", model.type = "spline")
PlotMiQC(pbmc3k, color.by = "miQC.keep")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/miQC_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Also, *RunMiQC* defaults to removing any cell with 75% or greater
posterior probability of being compromised, but if we want to be more or
less stringent, we can alter the *posterior.cutoff* parameter, like so:

``` r
pbmc3k <- RunMiQC(pbmc3k, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.9, 
    model.slot = "flexmix_model")
PlotMiQC(pbmc3k, color.by = "miQC.keep")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/miQC_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Note that when performing miQC multiple times on different samples for
the same experiment, it’s recommended to select the same
*posterior\_cutoff* for all, to give consistency in addition to the
flexibility of sample-specific models.

## When not to use miQC

The miQC model is based on the assumption that there are a non-trivial
number of compromised cells in the dataset, which is not true in all
datasets. We recommend using *FeatureScatter* on a dataset before
running miQC to see if the two-distribution model is appropriate. Look
for the distinctive triangular shape where cells have a wide variety of
mitochondrial percentages at lower gene counts and taper off to lower
mitochondrial percentage at higher gene counts.

For example of a dataset where there’s not a significant number of
compromised cells, so the two-distribution assumption is not met, we
simulate an extreme case using the `"pbmc3k"` dataset here.

``` r
set.seed(2021)
pbmc3k_extreme <- pbmc3k
simulated_percent_mt <- rnorm(mean = 2.5, sd = 0.2, n = ncol(pbmc3k_extreme))
pbmc3k_extreme$percent.mt <- ifelse(pbmc3k_extreme$nFeature_RNA > 400, simulated_percent_mt, pbmc3k_extreme$percent.mt)
simulated_percent_mt_2 <- runif(min = 0, max = 60, n = ncol(pbmc3k_extreme))
pbmc3k_extreme$percent.mt <- ifelse(pbmc3k_extreme$nFeature_RNA < 400, simulated_percent_mt_2, pbmc3k_extreme$percent.mt)
FeatureScatter(pbmc3k_extreme, feature1 = "nFeature_RNA", feature2 = "percent.mt")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/miQC_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

The *RunMiQC* function will throw a warning if only one distribution is
found. In these cases, we recommend using other filtering methods, such
as a cutoff on mitochondrial percentage or percentile using the
`"backup.option"` parameter to one of `"c("percentile", "percent",
"pass", "halt")`.

``` r
pbmc3k_extreme <- RunMiQC(pbmc3k_extreme, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", 
    posterior.cutoff = 0.9, model.slot = "flexmix_model", backup.option = "percentile", backup.percentile = 0.95)
```

    ## Warning in RunMiQC(pbmc3k_extreme, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", : flexmix returned only 1
    ## cluster

    ## defaulting to backup.percentile for filtering

    ## Warning: Adding a command log without an assay associated with it

``` r
FeatureScatter(pbmc3k_extreme, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "miQC.keep")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/miQC_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Session Information

    ## R version 4.0.4 (2021-02-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] pbmc3k.SeuratData_3.1.4 flexmix_2.3-17          lattice_0.20-41         SeuratWrappers_0.3.0   
    ## [5] SeuratData_0.2.1        SeuratObject_4.0.1      Seurat_4.0.1           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15            colorspace_2.0-2      deldir_0.2-10         modeltools_0.2-23     ellipsis_0.3.2       
    ##   [6] ggridges_0.5.3        rprojroot_2.0.2       spatstat.data_2.1-0   farver_2.1.0          leiden_0.3.7         
    ##  [11] listenv_0.8.0         remotes_2.3.0         ggrepel_0.9.1         fansi_0.5.0           R.methodsS3_1.8.1    
    ##  [16] codetools_0.2-18      splines_4.0.4         knitr_1.33            polyclip_1.10-0       jsonlite_1.7.2       
    ##  [21] ica_1.0-2             cluster_2.1.0         R.oo_1.24.0           png_0.1-7             uwot_0.1.10          
    ##  [26] shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0 BiocManager_1.30.15   compiler_4.0.4       
    ##  [31] httr_1.4.2            Matrix_1.3-3          fastmap_1.1.0         lazyeval_0.2.2        cli_3.0.1            
    ##  [36] later_1.2.0           formatR_1.9           htmltools_0.5.1.1     prettyunits_1.1.1     tools_4.0.4          
    ##  [41] rsvd_1.0.5            igraph_1.2.6          gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
    ##  [46] reshape2_1.4.4        dplyr_1.0.6           rappdirs_0.3.3        Rcpp_1.0.6            scattermore_0.7      
    ##  [51] jquerylib_0.1.4       vctrs_0.3.8           nlme_3.1-152          lmtest_0.9-38         xfun_0.23            
    ##  [56] stringr_1.4.0         globals_0.14.0        ps_1.6.0              mime_0.10             miniUI_0.1.1.1       
    ##  [61] lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2         future_1.21.0         MASS_7.3-53          
    ##  [66] zoo_1.8-9             scales_1.1.1          spatstat.core_2.1-2   promises_1.2.0.1      spatstat.utils_2.1-0 
    ##  [71] parallel_4.0.4        RColorBrewer_1.1-2    yaml_2.2.1            curl_4.3.1            reticulate_1.20      
    ##  [76] pbapply_1.4-3         gridExtra_2.3         ggplot2_3.3.5         sass_0.4.0            rpart_4.1-15         
    ##  [81] stringi_1.6.2         highr_0.9             pkgbuild_1.2.0        rlang_0.4.11          pkgconfig_2.0.3      
    ##  [86] matrixStats_0.59.0    evaluate_0.14         tensor_1.5            ROCR_1.0-11           purrr_0.3.4          
    ##  [91] labeling_0.4.2        patchwork_1.1.1       htmlwidgets_1.5.3     cowplot_1.1.1         processx_3.5.2       
    ##  [96] tidyselect_1.1.1      parallelly_1.25.0     RcppAnnoy_0.0.18      plyr_1.8.6            magrittr_2.0.1       
    ## [101] R6_2.5.0              generics_0.1.0        mgcv_1.8-33           pillar_1.6.1          withr_2.4.2          
    ## [106] fitdistrplus_1.1-3    nnet_7.3-15           abind_1.4-5           survival_3.2-7        tibble_3.1.2         
    ## [111] future.apply_1.7.0    crayon_1.4.1          KernSmooth_2.23-18    utf8_1.2.1            spatstat.geom_2.1-0  
    ## [116] plotly_4.9.3          rmarkdown_2.8         grid_4.0.4            data.table_1.14.0     callr_3.7.0          
    ## [121] digest_0.6.27         xtable_1.8-4          tidyr_1.1.3           httpuv_1.6.1          R.utils_2.10.1       
    ## [126] stats4_4.0.4          munsell_0.5.0         viridisLite_0.4.0     bslib_0.2.5.1
