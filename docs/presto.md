Fast Differential Expression with Presto
================
Compiled: October 07, 2020

This vignette demonstrates the use of the Presto package in Seurat.
Commands and parameters are based off of the [Presto
tutorial](http://htmlpreview.github.io/?https://github.com/immunogenomics/presto/blob/master/docs/getting-started.html).
If you use Presto in your work, please cite:

> *Presto scales Wilcoxon and auROC analyses to millions of
> observations*
> 
> Ilya Korsunsky, Aparna Nathan, Nghia Millard, Soumya Raychaudhuri
> 
> bioRxiv, 2019.
> 
> Pre-print: <https://www.biorxiv.org/content/10.1101/653253v1.full.pdf>
> 
> GitHub: <https://github.com/immunogenomics/presto>

Prerequisites to install:

  - [Seurat](https://satijalab.org/seurat/install)
  - [Presto](https://github.com/immunogenomics/presto)
  - [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
  - [SeuratData](https://github.com/satijalab/seurat-data)

<!-- end list -->

``` r
library(presto)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
```

### Differential Expression Testing for PBMC scRNA-seq Data

To learn more about this dataset, type `?pbmc3k`

``` r
InstallData("pbmc3k")
data("pbmc3k")
pbmc3k <- NormalizeData(pbmc3k)
Idents(pbmc3k) <- "seurat_annotations"

diffexp.B.Mono <- RunPresto(pbmc3k, "CD14+ Mono", "B")
head(diffexp.B.Mono, 10)
```

    ##                p_val avg_logFC pct.1 pct.2     p_val_adj
    ## CD79A  1.660326e-143 -2.989854 0.042 0.936 2.276972e-139
    ## TYROBP 3.516407e-138  3.512505 0.994 0.102 4.822401e-134
    ## S100A9 7.003189e-137  4.293303 0.996 0.134 9.604174e-133
    ## CST3   1.498348e-135  3.344758 0.992 0.174 2.054834e-131
    ## S100A4 8.872946e-135  2.854897 1.000 0.360 1.216836e-130
    ## LYZ    2.720838e-134  3.788514 1.000 0.422 3.731357e-130
    ## S100A8 3.115452e-133  4.039777 0.975 0.076 4.272530e-129
    ## CD79B  8.317731e-133 -2.667534 0.083 0.916 1.140694e-128
    ## S100A6 5.156920e-132  2.541609 0.996 0.352 7.072201e-128
    ## LGALS1 1.427548e-131  3.002493 0.979 0.131 1.957739e-127

``` r
diffexp.all <- RunPrestoAll(pbmc3k)
head(diffexp.all[diffexp.all$cluster == "B", ], 10)
```

    ##                     p_val avg_logFC pct.1 pct.2     p_val_adj cluster      gene
    ## CD79A.3      0.000000e+00  2.933865 0.936 0.044  0.000000e+00       B     CD79A
    ## MS4A1.3      0.000000e+00  2.290577 0.855 0.055  0.000000e+00       B     MS4A1
    ## LINC00926.1 2.998236e-274  1.956493 0.564 0.010 4.111781e-270       B LINC00926
    ## CD79B.3     1.126919e-273  2.381160 0.916 0.144 1.545457e-269       B     CD79B
    ## TCL1A.3     1.962618e-272  2.463556 0.622 0.023 2.691534e-268       B     TCL1A
    ## HLA-DQA1.2  3.017803e-267  2.104207 0.890 0.119 4.138616e-263       B  HLA-DQA1
    ## VPREB3      2.131575e-238  1.667466 0.488 0.008 2.923242e-234       B    VPREB3
    ## HLA-DQB1.2  2.076231e-230  2.112052 0.863 0.148 2.847343e-226       B  HLA-DQB1
    ## CD74.2      1.000691e-184  2.010688 1.000 0.819 1.372347e-180       B      CD74
    ## HLA-DRA.3   1.813356e-184  1.914531 1.000 0.492 2.486837e-180       B   HLA-DRA
