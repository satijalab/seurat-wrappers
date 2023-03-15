Here we illustrate the use of the
[fastTopics](https://github.com/stephenslab/fastTopics) Seurat wrapper
to analyze a Seurat data set. This vignette is only intended to
introduce the basic fastTopics interface for Seurat objects—for
background and guidance on analysis of single-cell RNA-seq data using a
topic model, please see the [fastTopics
vignettes](https://stephenslab.github.io/fastTopics/articles).

If you find the **fastTopics** package useful for your work, please
cite:

K. K. Dey, C. Joyce Hsiao and M. Stephens (2017). [Visualizing the
structure of RNA-seq expression data using grade of membership
models.](https://doi.org/10.1371/journal.pgen.1006599) PLoS Genetics 13,
e1006599.

P. Carbonetto, A. Sarkar, Z. Wang and M. Stephens (2021). [Non-negative
matrix factorization algorithms greatly improve topic model
fits.](https://arxiv.org/abs/2105.13440) arXiv 2105.13440

If you used the `de_analysis` function in fastTopics, please cite:

P. Carbonetto, K. Luo, A. Sarkar, A. Hung, K. Tayeb, S. Pott and M.
Stephens (2023). [Interpreting structure in sequence count data with
differential expression analysis allowing for grades of
membership.](https://doi.org/10.1101/2023.03.03.531029) bioRxiv
<doi:10.1101/2023.03.03.531029>

We begin by loading the packages used to perform the analysis.

    library(Seurat)
    library(SeuratData)
    library(SeuratWrappers)
    library(fastTopics)
    library(cowplot)

We set the seed so that the results can be reproduced.

    set.seed(1)

Load—and, if necessary, install—the PBMC 3k data set containing
transcription profiles for 2,700 cells.

    InstallData("pbmc3k")
    data(pbmc3k)
    dim(GetAssayData(pbmc3k))
    # [1] 13714  2700

Fit the multinomial topic model to the raw UMI counts—*no pre-processing
or pre-selection of genes is needed.* Note that it may take several
minutes to complete this model fitting step.

    pbmc3k <- FitTopicModel(pbmc3k,k = 6)

To fit a topic model, we must specify *K*, the number of topics. Here,
we have chosen *K* = 6 topics. In most settings, a good choice of *K*
will not be known in advance, so you will you want to explore the
results from topic models at different settings of *K*.

This plot shows the cells projected onto the top two principal
components (PCs) of the topic mixture proportions.

    Idents(pbmc3k) <- pbmc3k$seurat_annotations
    DimPlot(pbmc3k,reduction = "pca_topics",pt.size = 1) +
      theme_cowplot(font_size = 10)

<img src="fasttopics_files/figure-markdown_strict/pca-1-1.png" style="display: block; margin: auto;" />

Compare this against the top two PCs of the transformed counts:

    pbmc3k <- FindVariableFeatures(pbmc3k)
    pbmc3k <- NormalizeData(pbmc3k)
    pbmc3k <- ScaleData(pbmc3k)
    pbmc3k <- RunPCA(pbmc3k)
    DimPlot(pbmc3k,reduction = "pca",pt.size = 1) +
      theme_cowplot(font_size = 10)

<img src="fasttopics_files/figure-markdown_strict/pca-2-1.png" style="display: block; margin: auto;" />

The fitted topic model—a “multinom\_topic\_model” object—is stored in
the “misc” slot:

    fit <- Misc(Reductions(pbmc3k,"multinom_topic_model"))

Once the fitted topic model is extracted, many functions from the
**fastTopics** package can be used for analysis and visualization. For
example, the Structure plot provides an evocative visual summary of the
estimated mixture proportions for each cell. Here, we have grouped the
cells by previously inferred labels.

    structure_plot(fit,grouping = Idents(pbmc3k),gap = 25)

<img src="fasttopics_files/figure-markdown_strict/structure-plot-1.png" style="display: block; margin: auto;" />

This is the version of R and the packages that were used to generate
these results.

    sessionInfo()
    # R version 4.2.0 (2022-04-22)
    # Platform: x86_64-pc-linux-gnu (64-bit)
    # Running under: Red Hat Enterprise Linux 8.4 (Ootpa)
    # 
    # Matrix products: default
    # BLAS/LAPACK: /software/openblas-0.3.13-el8-x86_64/lib/libopenblas_skylakexp-r0.3.13.so
    # 
    # locale:
    #  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    #  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    #  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    # [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    # 
    # attached base packages:
    # [1] stats     graphics  grDevices utils     datasets  methods   base     
    # 
    # other attached packages:
    # [1] rmarkdown_2.14          SeuratWrappers_0.3.1    pbmc3k.SeuratData_3.1.4
    # [4] cowplot_1.1.1           fastTopics_0.6-150      SeuratData_0.2.2       
    # [7] SeuratObject_4.1.3      Seurat_4.3.0           
    # 
    # loaded via a namespace (and not attached):
    #   [1] utf8_1.2.2             spatstat.explore_3.0-6 reticulate_1.24       
    #   [4] R.utils_2.11.0         tidyselect_1.1.2       htmlwidgets_1.5.4     
    #   [7] grid_4.2.0             Rtsne_0.16             devtools_2.4.3        
    #  [10] munsell_0.5.0          codetools_0.2-18       ica_1.0-2             
    #  [13] future_1.25.0          miniUI_0.1.1.1         withr_2.5.0           
    #  [16] spatstat.random_3.1-3  colorspace_2.0-3       progressr_0.10.0      
    #  [19] highr_0.9              knitr_1.39             rstudioapi_0.13       
    #  [22] ROCR_1.0-11            tensor_1.5             listenv_0.8.0         
    #  [25] labeling_0.4.2         mixsqp_0.3-48          polyclip_1.10-0       
    #  [28] MCMCpack_1.6-3         farver_2.1.0           rprojroot_2.0.3       
    #  [31] coda_0.19-4            parallelly_1.31.1      vctrs_0.4.1           
    #  [34] generics_0.1.2         xfun_0.30              R6_2.5.1              
    #  [37] rsvd_1.0.5             invgamma_1.1           spatstat.utils_3.0-1  
    #  [40] cachem_1.0.6           assertthat_0.2.1       promises_1.2.0.1      
    #  [43] scales_1.2.0           gtable_0.3.0           globals_0.14.0        
    #  [46] processx_3.5.3         goftest_1.2-3          mcmc_0.9-7            
    #  [49] rlang_1.0.2            MatrixModels_0.5-0     splines_4.2.0         
    #  [52] lazyeval_0.2.2         spatstat.geom_3.0-6    BiocManager_1.30.20   
    #  [55] yaml_2.3.5             reshape2_1.4.4         abind_1.4-5           
    #  [58] httpuv_1.6.5           tools_4.2.0            usethis_2.1.5         
    #  [61] ggplot2_3.3.6          ellipsis_0.3.2         jquerylib_0.1.4       
    #  [64] RColorBrewer_1.1-3     sessioninfo_1.2.2      ggridges_0.5.3        
    #  [67] Rcpp_1.0.9             plyr_1.8.7             progress_1.2.2        
    #  [70] purrr_0.3.4            ps_1.7.0               prettyunits_1.1.1     
    #  [73] deldir_1.0-6           pbapply_1.5-0          ashr_2.2-54           
    #  [76] zoo_1.8-10             ggrepel_0.9.1          cluster_2.1.3         
    #  [79] fs_1.5.2               magrittr_2.0.3         data.table_1.14.4     
    #  [82] scattermore_0.8        SparseM_1.81           lmtest_0.9-40         
    #  [85] RANN_2.6.1             truncnorm_1.0-8        SQUAREM_2021.1        
    #  [88] fitdistrplus_1.1-8     matrixStats_0.62.0     pkgload_1.2.4         
    #  [91] hms_1.1.1              patchwork_1.1.1        mime_0.12             
    #  [94] evaluate_0.15          xtable_1.8-4           gridExtra_2.3         
    #  [97] testthat_3.1.4         compiler_4.2.0         tibble_3.1.7          
    # [100] KernSmooth_2.23-20     crayon_1.5.1           R.oo_1.24.0           
    # [103] htmltools_0.5.2        later_1.3.0            tidyr_1.2.0           
    # [106] RcppParallel_5.1.5     DBI_1.1.2              MASS_7.3-56           
    # [109] rappdirs_0.3.3         Matrix_1.5-3           brio_1.1.3            
    # [112] cli_3.3.0              quadprog_1.5-8         R.methodsS3_1.8.1     
    # [115] parallel_4.2.0         igraph_1.3.1           pkgconfig_2.0.3       
    # [118] sp_1.6-0               plotly_4.10.0          spatstat.sparse_3.0-0 
    # [121] xml2_1.3.3             roxygen2_7.1.2         bslib_0.3.1           
    # [124] stringr_1.4.0          callr_3.7.0            digest_0.6.29         
    # [127] sctransform_0.3.5      RcppAnnoy_0.0.19       spatstat.data_3.0-0   
    # [130] leiden_0.3.10          uwot_0.1.14            shiny_1.7.1           
    # [133] quantreg_5.93          lifecycle_1.0.1        nlme_3.1-157          
    # [136] jsonlite_1.8.0         desc_1.4.1             viridisLite_0.4.0     
    # [139] fansi_1.0.3            pillar_1.7.0           lattice_0.20-45       
    # [142] fastmap_1.1.0          httr_1.4.2             pkgbuild_1.3.1        
    # [145] survival_3.3-1         glue_1.6.2             remotes_2.4.2         
    # [148] png_0.1-7              stringi_1.7.6          sass_0.4.1            
    # [151] memoise_2.0.1          dplyr_1.0.9            irlba_2.3.5           
    # [154] future.apply_1.9.0
