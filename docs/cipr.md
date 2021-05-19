Using CIPR with human PBMC data
================
Atakan Ekiz
19 May, 2021

This vignette demonstrates how to run CIPR on Seurat objects. If you use
CIPR, please cite:

> *CIPR: a web-based R/shiny app and R package to annotate cell clusters
> in single cell RNA sequencing experiments*
> 
> H. Atakan Ekiz, Christopher J. Conley, W. Zac Stephens & Ryan M.
> O’Connell
> 
> BMC Bioinformatics, 2020.
> 
> doi:
> [10.1186/s12859-020-3538-2](https://doi.org/10.1186/s12859-020-3538-2)
> 
> Github: <https://github.com/atakanekiz/CIPR-Package>

# Summary

This vignette describes how to use CIPR package with 3k PBMC data freely
available from 10X genomics. Here, we recycle the code described in
[Seurat’s guided clustering
tutorial](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html) to
help users perform analyses from scratch. Using this dataset we will
demonstrate the capabilities of CIPR to annotate single cell clusters in
single cell RNAseq (scRNAseq) experiments. For further information about
other clustering methods, please see Seurat’s comprehensive
[website](https://satijalab.org/seurat/)

# Install CIPR

``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Use this option if you want to build vignettes during installation This can take a long time
# due to the installation of suggested packages.
remotes::install_github("atakanekiz/CIPR-Package", build_vignettes = TRUE)

# Use this if you would like to install the package without vignettes
# remotes::install_github('atakanekiz/CIPR-Package')
```

# Seurat pipeline

## Setup Seurat object

``` r
library(dplyr)
library(Seurat)
library(SeuratData)
library(CIPR)
```

``` r
# Load data
InstallData("pbmc3k")
pbmc <- pbmc3k
```

## Pre-processing

The steps below encompass the standard pre-processing workflow for
scRNA-seq data in Seurat. These represent the selection and filtration
of cells based on QC metrics, data normalization and scaling, and the
detection of highly variable features.

``` r
# Calculate mitochondrial gene representation (indicative of low quality cells)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter out genes with feature counts outside of 200-2500 range, and >5% mt genes
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

## Normalizing data

``` r
pbmc <- NormalizeData(pbmc)
```

## Variable gene detection and scaling

``` r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

## Perform PCA

``` r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

``` r
ElbowPlot(pbmc)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Cluster cells

``` r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

## Run non-linear dimensionality reduction (tSNE)

``` r
pbmc <- RunTSNE(pbmc, dims = 1:10)
pbmc$unnamed_clusters <- Idents(pbmc)
```

``` r
# saveRDS(pbmc, 'pbmc.rds')
```

## Find differentially expressed genes

**This is the step where we generate the input for CIPR’s log fold
change (logFC) comparison methods.**

## Calculate average gene expression per cluster

**This is the step where we generate the input for CIPR’s all-genes
correlation methods.**

``` r
avgexp <- AverageExpression(pbmc)
avgexp <- as.data.frame(x = avgexp$RNA)
avgexp$gene <- rownames(avgexp)
```

## Visualize Seurat pbject

``` r
DimPlot(pbmc)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

# CIPR analysis

The user can select one of the 7 provided reference data sets:

| Reference                                 | `reference` argument |
| ----------------------------------------- | -------------------- |
| Immunological Genome Project (ImmGen)     | “immgen”             |
| Presorted cell RNAseq (various tissues)   | “mmrnaseq”           |
| Blueprint/ENCODE                          | “blueprint”          |
| Human Primary Cell Atlas                  | “hpca”               |
| Database of Immune Cell Expression (DICE) | “dice”               |
| Hematopoietic differentiation             | “hema”               |
| Presorted cell RNAseq (PBMC)              | “hsrnaseq”           |
| User-provided custom reference            | “custom”             |

## Standard logFC comparison method

In this method CIPR accepts `allmarkers` data frame created above and
performs the following analytical steps:

  - It calculates a vector of logFC values for each reference sample
    (i.e. individual columns of the reference data frame) by comparing
    log-normalized expression value of a gene (i.e. rows of the
    reference data frame) to the average gene expression across the
    entire reference dataset.
  - It then scores unknown cluster logFC differential gene expression
    data against this reference logFC values to create a vector of
    identity scores
  - User can select one of three methods:
      - LogFC dot product (sum of all logFC x logFC values among
        matching genes). This is the recommended method in CIPR.
      - LogFC Spearman’s correlation (rank correlation of logFC values)
      - LogFC Pearson’s correlation (linear correlation of logFC values)

### Plot all identity scores per cluster-reference cell pairs

The code below performs analysis using sorted human PBMC RNAseq data as
reference, and plots

CIPR results can be summarized for each cluster in scatter plots.

``` r
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
     # axis.text.x=element_text(color="red") # arguments to pass to ggplot2::theme() to change plotting parameters
     )
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### Plot identity scores for a select cluster

`ind_clu_plots` object is created in the global environment to help
users can visualize results for a desired cluster and manipulate
graphing parameters. ggplot2 functions can be iteratively added to
individual plots to create annotations etc.

``` r
library(ggplot2)
ind_clu_plots$cluster6 + theme(axis.text.y = element_text(color = "red"), axis.text.x = element_text(color = "blue")) + 
    labs(fill = "Reference") + ggtitle("Figure S4a. Automated cluster annotation results are shown for cluster 6") + 
    annotate("text", label = "2 sd range", x = 10, y = 700, size = 8, color = "steelblue") + annotate("text", 
    label = "1 sd range", x = 10, y = 200, size = 8, color = "orange2") + geom_rect(aes(xmin = 94, 
    xmax = 99, ymin = 1000, ymax = 1300), fill = NA, size = 3, color = "red")
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

### Plot top scoring refernce subsets for each cluster

``` r
CIPR(input_dat = allmarkers, comp_method = "logfc_dot_product", reference = "hsrnaseq", plot_ind = F, 
    plot_top = T, global_results_obj = T, global_plot_obj = T)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Tabulate CIPR results

CIPR results (both top 5 scoring reference types per cluster and the
entire analysis) are saved as global objects (`CIPR_top_results` and
`CIPR_all_results` respectively) to allow users to explore the outputs
and generate specific plots and tables.

``` r
head(CIPR_top_results)
```

    ## # A tibble: 6 x 9
    ## # Groups:   cluster [2]
    ##   cluster reference_cell_t… reference_id  long_name   description identity_score
    ##   <fct>   <fct>             <chr>         <chr>       <chr>                <dbl>
    ## 1 0       CD8+ T cell       G4YW_CD8_nai… Naive CD8 … N/A                   838.
    ## 2 0       CD8+ T cell       DZQV_CD8_nai… Naive CD8 … N/A                   833.
    ## 3 0       CD8+ T cell       925L_CD8_nai… Naive CD8 … N/A                   779.
    ## 4 0       CD8+ T cell       9JD4_CD8_nai… Naive CD8 … N/A                   751.
    ## 5 0       CD4+ T cell       9JD4_CD4_nai… Naive CD4 … N/A                   743.
    ## 6 1       Monocyte          G4YW_C_mono   Classical … N/A                  2031.
    ## # … with 3 more variables: index <int>, z_score <dbl>,
    ## #   percent_pos_correlation <dbl>

``` r
head(CIPR_all_results)
```

    ##        reference_id identity_score reference_cell_type
    ## 1      DZQV_B_naive      -506.4224              B cell
    ## 2        DZQV_B_NSM      -414.3927              B cell
    ## 3         DZQV_B_Ex      -438.5500              B cell
    ## 4         DZQV_B_SM      -441.4376              B cell
    ## 5 DZQV_Plasmablasts       226.2113              B cell
    ## 6      925L_B_naive      -128.9296              B cell
    ##                     long_name description cluster    z_score
    ## 1               Naive B cells         N/A       0 -1.0021725
    ## 2 Non-switched memory B cells         N/A       0 -0.8200527
    ## 3           Exhausted B cells         N/A       0 -0.8678580
    ## 4     Switched memory B cells         N/A       0 -0.8735724
    ## 5                Plasmablasts         N/A       0  0.4476555
    ## 6               Naive B cells         N/A       0 -0.2551421
    ##   percent_pos_correlation
    ## 1                42.16336
    ## 2                43.37748
    ## 3                43.92936
    ## 4                41.39073
    ## 5                64.56954
    ## 6                61.92053

## Standard all-genes correlation method

CIPR also implements a simple correlation approach in which overall
correlation in gene expression is calculated for the pairs of unknown
clusters and the reference samples (regardless of the differential
expression status of the gene). This approach is conceptually similar to
some other automated identity prediction pipelines such as
[SingleR](https://www.ncbi.nlm.nih.gov/pubmed/30643263) and
[scMCA](https://www.ncbi.nlm.nih.gov/pubmed/30758821).

  - **Spearman’s correlation:** It calculates correlation based on
    ranked gene expression. It can be suitable for comparing
    experimental and reference data which were obtained using different
    technologies.
  - **Pearson’s correlation:** It calculates linear correlations. This
    can be useful when the user would like to provide a custom reference
    dataset to CIPR.

### Plot all identity scores per cluster-reference cell pairs

The code below performs analysis using sorted human PBMC RNAseq data as
reference, and plots

CIPR results can be summarized for each cluster in scatter plots.

``` r
CIPR(input_dat = avgexp, comp_method = "all_genes_spearman", reference = "hsrnaseq", plot_ind = T, 
    plot_top = F, global_results_obj = T, global_plot_obj = T)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### Plot top scoring refernce subsets for each cluster

``` r
CIPR(input_dat = avgexp, comp_method = "all_genes_spearman", reference = "hsrnaseq", plot_ind = F, 
    plot_top = T, global_results_obj = T, global_plot_obj = T)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

### Tabulate CIPR results

CIPR results (both top 5 scoring reference types per cluster and the
entire analysis) are saved as global objects (`CIPR_top_results` and
`CIPR_all_results` respectively) to allow users to explore the outputs
and generate specific plots and tables.

``` r
head(CIPR_top_results)
```

    ## # A tibble: 6 x 8
    ## # Groups:   cluster [2]
    ##   cluster reference_cell_t… reference_id long_name    description identity_score
    ##   <fct>   <fct>             <chr>        <chr>        <chr>                <dbl>
    ## 1 0       CD4+ T cell       DZQV_CD4_na… Naive CD4 T… N/A                  0.797
    ## 2 0       CD4+ T cell       925L_TFH     Follicular … N/A                  0.793
    ## 3 0       CD4+ T cell       G4YW_Th1     Th1 cells    N/A                  0.788
    ## 4 0       CD4+ T cell       G4YW_Treg    T regulator… N/A                  0.786
    ## 5 0       CD4+ T cell       DZQV_Th17    Th17 cells   N/A                  0.780
    ## 6 1       Monocyte          G4YW_I_mono  Intermediat… N/A                  0.784
    ## # … with 2 more variables: index <int>, z_score <dbl>

``` r
head(CIPR_all_results)
```

    ##        reference_id identity_score reference_cell_type
    ## 1      DZQV_B_naive      0.6503197              B cell
    ## 2        DZQV_B_NSM      0.6480390              B cell
    ## 3         DZQV_B_Ex      0.6488979              B cell
    ## 4         DZQV_B_SM      0.6961983              B cell
    ## 5 DZQV_Plasmablasts      0.6816080              B cell
    ## 6      925L_B_naive      0.6421836              B cell
    ##                     long_name description cluster      z_score
    ## 1               Naive B cells         N/A       0 -0.636484838
    ## 2 Non-switched memory B cells         N/A       0 -0.668063024
    ## 3           Exhausted B cells         N/A       0 -0.656170534
    ## 4     Switched memory B cells         N/A       0 -0.001236968
    ## 5                Plasmablasts         N/A       0 -0.203258085
    ## 6               Naive B cells         N/A       0 -0.749138568

## Limiting analysis to the select subsets of reference data

Sometimes excluding irrelevant reference cell types from the analysis
can be helpful. Especially when the logFC comparison methods are
utilized, removing irrelevant subsets may improve discrimination of
closely related subsets, since the reference logFC values will be
calculated after subsetting the data frame. Filtering out reference
subsets should not impact results of the all-genes correlation methods,
but it can make the graphical outputs easier to look at

3k PBMC dataset may not be the best example to demonstrate benefits of
reference dataset subsetting, but the code below serves as an example
for this functionality.

``` r
CIPR(input_dat = allmarkers, comp_method = "logfc_dot_product", reference = "hsrnaseq", plot_ind = T, 
    plot_top = F, global_results_obj = T, global_plot_obj = T, select_ref_subsets = c("CD4+ T cell", 
        "CD8+ T cell", "Monocyte", "NK cell"))
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## Filtering out lowly variable genes

Genes that have a low expression variance across the reference data
frame has weaker discriminatory potential. Thus, excluding these genes
from the analysis can reduce the noise and improve the prediction
scores, especially when using all-genes correlation based methods.

We implemented a variance filtering parameter, `keep_top_var`, which
allows users to keep top Nth% variable reference genes in the analysis.
For instance, by setting this argument to 10, CIPR can be instructed to
use only the top 10% highly variable genes in identity score
calculations. In our experience *(Ekiz HA, BMC Bioinformatics, in
revision)* limiting the analysis to highly variable genes does not
significantly impact the identity scores of the top-scoring reference
cell subsets, but it reduces the identity scores of
intermediate/low-scoring reference cells leading to an improvement of
z-scores. The “best” value for this parameter remains to be determined
by the user in individual studies.

``` r
CIPR(input_dat = avgexp, comp_method = "all_genes_spearman", reference = "hsrnaseq", plot_ind = T, 
    plot_top = F, global_results_obj = T, global_plot_obj = T, keep_top_var = 10)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/cipr_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

<details>

<summary>**Session Info**</summary>

``` r
sessionInfo()
```

    ## R version 4.0.4 (2021-02-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_3.3.3           pbmc3k.SeuratData_3.1.4 CIPR_0.1.0             
    ## [4] SeuratData_0.2.1        SeuratObject_4.0.1      Seurat_4.0.1           
    ## [7] dplyr_1.0.6            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1          backports_1.2.1       plyr_1.8.6           
    ##   [4] igraph_1.2.6          lazyeval_0.2.2        splines_4.0.4        
    ##   [7] listenv_0.8.0         scattermore_0.7       digest_0.6.27        
    ##  [10] htmltools_0.5.1.1     fansi_0.4.2           magrittr_2.0.1       
    ##  [13] tensor_1.5            cluster_2.1.0         ROCR_1.0-11          
    ##  [16] openxlsx_4.2.3        limma_3.46.0          remotes_2.3.0        
    ##  [19] globals_0.14.0        matrixStats_0.58.0    spatstat.sparse_2.0-0
    ##  [22] prettyunits_1.1.1     colorspace_2.0-1      rappdirs_0.3.3       
    ##  [25] ggrepel_0.9.1         haven_2.4.1           xfun_0.23            
    ##  [28] callr_3.7.0           crayon_1.4.1          jsonlite_1.7.2       
    ##  [31] spatstat.data_2.1-0   survival_3.2-7        zoo_1.8-9            
    ##  [34] glue_1.4.2            polyclip_1.10-0       gtable_0.3.0         
    ##  [37] leiden_0.3.7          car_3.0-10            pkgbuild_1.2.0       
    ##  [40] future.apply_1.7.0    abind_1.4-5           scales_1.1.1         
    ##  [43] rstatix_0.7.0         miniUI_0.1.1.1        Rcpp_1.0.6           
    ##  [46] viridisLite_0.4.0     xtable_1.8-4          reticulate_1.20      
    ##  [49] spatstat.core_2.1-2   foreign_0.8-81        htmlwidgets_1.5.3    
    ##  [52] httr_1.4.2            RColorBrewer_1.1-2    ellipsis_0.3.2       
    ##  [55] ica_1.0-2             pkgconfig_2.0.3       farver_2.1.0         
    ##  [58] uwot_0.1.10           deldir_0.2-10         utf8_1.2.1           
    ##  [61] tidyselect_1.1.1      labeling_0.4.2        rlang_0.4.11         
    ##  [64] reshape2_1.4.4        later_1.2.0           cellranger_1.1.0     
    ##  [67] munsell_0.5.0         tools_4.0.4           cli_2.5.0            
    ##  [70] generics_0.1.0        broom_0.7.6           ggridges_0.5.3       
    ##  [73] evaluate_0.14         stringr_1.4.0         fastmap_1.1.0        
    ##  [76] yaml_2.2.1            goftest_1.2-2         processx_3.5.2       
    ##  [79] knitr_1.33            fitdistrplus_1.1-3    zip_2.1.1            
    ##  [82] purrr_0.3.4           RANN_2.6.1            pbapply_1.4-3        
    ##  [85] future_1.21.0         nlme_3.1-152          mime_0.10            
    ##  [88] formatR_1.9           compiler_4.0.4        rstudioapi_0.13      
    ##  [91] plotly_4.9.3          curl_4.3.1            png_0.1-7            
    ##  [94] ggsignif_0.6.1        spatstat.utils_2.1-0  tibble_3.1.2         
    ##  [97] stringi_1.6.2         highr_0.9             ps_1.6.0             
    ## [100] forcats_0.5.1         lattice_0.20-41       Matrix_1.3-3         
    ## [103] vctrs_0.3.8           pillar_1.6.1          lifecycle_1.0.0      
    ## [106] spatstat.geom_2.1-0   lmtest_0.9-38         RcppAnnoy_0.0.18     
    ## [109] data.table_1.14.0     cowplot_1.1.1         irlba_2.3.3          
    ## [112] httpuv_1.6.1          patchwork_1.1.1       R6_2.5.0             
    ## [115] promises_1.2.0.1      KernSmooth_2.23-18    gridExtra_2.3        
    ## [118] rio_0.5.26            parallelly_1.25.0     codetools_0.2-18     
    ## [121] MASS_7.3-53           gtools_3.8.2          rprojroot_2.0.2      
    ## [124] withr_2.4.2           sctransform_0.3.2     hms_1.1.0            
    ## [127] mgcv_1.8-33           parallel_4.0.4        grid_4.0.4           
    ## [130] rpart_4.1-15          tidyr_1.1.3           rmarkdown_2.8        
    ## [133] carData_3.0-4         Rtsne_0.15            ggpubr_0.4.0         
    ## [136] shiny_1.6.0

</details>
