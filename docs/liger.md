Integrating Seurat objects using LIGER
================
Compiled: May 25, 2021

  - [](#section)
      - [Systematic comparative analysis of human
        PBMC](#systematic-comparative-analysis-of-human-pbmc)
      - [Interferon-stimulated and control
        PBMC](#interferon-stimulated-and-control-pbmc)
      - [Eight human pancreatic islet
        datasets](#eight-human-pancreatic-islet-datasets)

NOTE: Please update your `liger` version to 0.5.0 or above before
following this tutorial.

This vigettte demonstrates how to run LIGER on Seurat objects.
Parameters and commands are based on the [LIGER
tutorial](http://htmlpreview.github.io/?https://github.com/MacoskoLab/liger/blob/master/vignettes/Integrating_multi_scRNA_data.html).
If you use LIGER, please cite:

> *Single-Cell Multi-omic Integration Compares and Contrasts Features of
> Brain Cell Identity*
> 
> Joshua Welch, Velina Kozareva, Ashley Ferreira, Charles Vanderburg,
> Carly Martin, Evan Z.Macosko
> 
> Cell, 2019.
> 
> doi:
> [10.1016/j.cell.2019.05.006](https://doi.org/10.1016/j.cell.2019.05.006)
> 
> GitHub: <https://github.com/MacoskoLab/liger>

Prerequisites to install:

  - [Seurat](https://satijalab.org/seurat/install)
  - [LIGER](https://github.com/MacoskoLab/liger)
  - [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
  - [SeuratData](https://github.com/satijalab/seurat-data)

<!-- end list -->

``` r
library(rliger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
```

In order to replicate LIGER’s multi-dataset functionality, we will use
the `split.by` parameter to preprocess the Seurat object on subsets of
the data belonging to each dataset separately. Also, as LIGER does not
center data when scaling, we will skip that step as well.

`RunQuantileNorm` produces joint clusters, but users can also optionally
perform Louvain community detection (`FindNeighbors` and `FindClusters`)
on the integrated latent space from iNMF.

## 

### Systematic comparative analysis of human PBMC

To learn more about this dataset, type `?pbmcsca`

``` r
InstallData("pbmcsca")
data("pbmcsca")
# Please update your `liger` version to 0.5.0 or above before following this tutorial
pbmcsca <- NormalizeData(pbmcsca)
pbmcsca <- FindVariableFeatures(pbmcsca)
pbmcsca <- ScaleData(pbmcsca, split.by = "Method", do.center = FALSE)
pbmcsca <- RunOptimizeALS(pbmcsca, k = 20, lambda = 5, split.by = "Method")
pbmcsca <- RunQuantileNorm(pbmcsca, split.by = "Method")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs
pbmcsca <- FindNeighbors(pbmcsca, reduction = "iNMF", dims = 1:20)
pbmcsca <- FindClusters(pbmcsca, resolution = 0.3)
# Dimensional reduction and plotting
pbmcsca <- RunUMAP(pbmcsca, dims = 1:ncol(pbmcsca[["iNMF"]]), reduction = "iNMF")
DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/liger_files/figure-gfm/pbmcsca-1.png)<!-- -->

### Interferon-stimulated and control PBMC

To learn more about this dataset, type `?ifnb`

``` r
InstallData("ifnb")
data("ifnb")
# Please update your `liger` version to 0.5.0 or above before following this tutorial.
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb, split.by = "stim", do.center = FALSE)
ifnb <- RunOptimizeALS(ifnb, k = 20, lambda = 5, split.by = "stim")
ifnb <- RunQuantileNorm(ifnb, split.by = "stim")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs
ifnb <- FindNeighbors(ifnb, reduction = "iNMF", dims = 1:20)
ifnb <- FindClusters(ifnb, resolution = 0.55)
# Dimensional reduction and plotting
ifnb <- RunUMAP(ifnb, dims = 1:ncol(ifnb[["iNMF"]]), reduction = "iNMF")
DimPlot(ifnb, group.by = c("stim", "ident", "seurat_annotations"), ncol = 3)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/liger_files/figure-gfm/ifnb-1.png)<!-- -->

### Eight human pancreatic islet datasets

To learn more about this dataset, type `?panc8`

``` r
InstallData("panc8")
data("panc8")
# Please update your `liger` version to 0.5.0 or above before following this tutorial.
panc8 <- NormalizeData(panc8)
panc8 <- FindVariableFeatures(panc8)
panc8 <- ScaleData(panc8, split.by = "replicate", do.center = FALSE)
panc8 <- RunOptimizeALS(panc8, k = 20, lambda = 5, split.by = "replicate")
panc8 <- RunQuantileNorm(panc8, split.by = "replicate")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs
panc8 <- FindNeighbors(panc8, reduction = "iNMF", dims = 1:20)
panc8 <- FindClusters(panc8, resolution = 0.4)
# Dimensional reduction and plotting
panc8 <- RunUMAP(panc8, dims = 1:ncol(panc8[["iNMF"]]), reduction = "iNMF")
DimPlot(panc8, group.by = c("replicate", "ident", "celltype"), ncol = 3)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/liger_files/figure-gfm/pancreas-1.png)<!-- -->
