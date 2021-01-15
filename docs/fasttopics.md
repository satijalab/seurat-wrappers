Add title here
================
Compiled: January 15, 2021

*Add introductory text here.*

*Add citation here.*

We begin by loading the packages used to perform the analysis.

``` r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(fastTopics)
```

We set the seed so that the results can be reproduced.

``` r
set.seed(1)
```

Load—and, if necessary, install—the PBMC 3k data
    set.

``` r
InstallData("pbmc3k")
```

    ## Warning: The following packages are already installed and will not be
    ## reinstalled: pbmc3k

``` r
data(pbmc3k)
```

This is the version of R and the packages that were used to generate
these results.

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fastTopics_0.4-23       SeuratWrappers_0.3.2    pbmc3k.SeuratData_3.1.4
    ## [4] SeuratData_0.2.1        Seurat_3.2.3           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15            colorspace_1.4-1      deldir_0.1-29        
    ##   [4] ggridges_0.5.2        spatstat.data_1.4-3   leiden_0.3.3         
    ##   [7] listenv_0.8.0         MatrixModels_0.4-1    remotes_2.1.0        
    ##  [10] ggrepel_0.9.0         fansi_0.4.0           codetools_0.2-16     
    ##  [13] splines_3.6.2         knitr_1.26            polyclip_1.10-0      
    ##  [16] zeallot_0.1.0         jsonlite_1.6          mcmc_0.9-6           
    ##  [19] ica_1.0-2             cluster_2.1.0         png_0.1-7            
    ##  [22] uwot_0.1.10           shiny_1.4.0           sctransform_0.3.2    
    ##  [25] BiocManager_1.30.10   compiler_3.6.2        httr_1.4.2           
    ##  [28] backports_1.1.5       assertthat_0.2.1      Matrix_1.2-18        
    ##  [31] fastmap_1.0.1         lazyeval_0.2.2        cli_2.0.0            
    ##  [34] later_1.0.0           htmltools_0.4.0       quantreg_5.54        
    ##  [37] prettyunits_1.1.1     tools_3.6.2           rsvd_1.0.2           
    ##  [40] igraph_1.2.5          coda_0.19-3           gtable_0.3.0         
    ##  [43] glue_1.3.1            RANN_2.6.1            reshape2_1.4.3       
    ##  [46] dplyr_0.8.3           rappdirs_0.3.1        Rcpp_1.0.5           
    ##  [49] spatstat_1.64-1       scattermore_0.7       vctrs_0.2.1          
    ##  [52] nlme_3.1-142          lmtest_0.9-38         xfun_0.11            
    ##  [55] stringr_1.4.0         globals_0.13.0        mime_0.8             
    ##  [58] miniUI_0.1.1.1        lifecycle_0.1.0       irlba_2.3.3          
    ##  [61] goftest_1.2-2         future_1.18.0         MASS_7.3-51.4        
    ##  [64] zoo_1.8-7             scales_1.1.0          hms_0.5.2            
    ##  [67] promises_1.1.0        spatstat.utils_1.17-0 parallel_3.6.2       
    ##  [70] SparseM_1.78          RColorBrewer_1.1-2    yaml_2.2.0           
    ##  [73] reticulate_1.16       pbapply_1.4-3         gridExtra_2.3        
    ##  [76] ggplot2_3.3.0         rpart_4.1-15          stringi_1.4.3        
    ##  [79] rlang_0.4.5           pkgconfig_2.0.3       matrixStats_0.56.0   
    ##  [82] evaluate_0.14         lattice_0.20-38       ROCR_1.0-11          
    ##  [85] purrr_0.3.3           tensor_1.5            patchwork_1.0.1      
    ##  [88] htmlwidgets_1.5.1     cowplot_1.0.0         tidyselect_0.2.5     
    ##  [91] RcppAnnoy_0.0.18      plyr_1.8.5            magrittr_1.5         
    ##  [94] R6_2.4.1              pillar_1.4.3          mgcv_1.8-31          
    ##  [97] fitdistrplus_1.1-1    survival_3.1-8        abind_1.4-5          
    ## [100] tibble_2.1.3          future.apply_1.6.0    crayon_1.3.4         
    ## [103] KernSmooth_2.23-16    plotly_4.9.2          rmarkdown_2.3        
    ## [106] progress_1.2.2        grid_3.6.2            data.table_1.12.8    
    ## [109] digest_0.6.23         xtable_1.8-4          tidyr_1.0.0          
    ## [112] httpuv_1.5.2          MCMCpack_1.4-5        RcppParallel_4.4.2   
    ## [115] munsell_0.5.0         viridisLite_0.3.0     quadprog_1.5-8
