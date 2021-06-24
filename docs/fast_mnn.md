Running fastMNN on Seurat Objects
================
Compiled: June 23, 2021

  - [](#section)
      - [Systematic comparative analysis of human
        PBMC](#systematic-comparative-analysis-of-human-pbmc)
      - [Interferon-stimulated and control
        PBMC](#interferon-stimulated-and-control-pbmc)
      - [Eight human pancreatic islet
        datasets](#eight-human-pancreatic-islet-datasets)

This vigettte demonstrates how to run fastMNN on Seurat objects.
Parameters and commands are based off of the [fastMNN help
page](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html). If you
use fastMNN, please cite:

> *Batch effects in single-cell RNA-sequencing data are corrected by
> matching mutual nearest neighbors*
> 
> Laleh Haghverdi, Aaron T L Lun, Michael D Morgan & John C Marioni
> 
> Nature Biotechnology, 2018
> 
> doi: [10.1038/nbt.4091](https://doi.org/10.1038/nbt.4091)
> 
> Bioconductor:
> <https://bioconductor.org/packages/release/bioc/html/batchelor.html>

Prerequisites to install:

  - [Seurat](https://satijalab.org/seurat/install)
  - [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
  - [SeuratData](https://github.com/satijalab/seurat-data)

<!-- end list -->

``` r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
```

## 

### Systematic comparative analysis of human PBMC

To learn more about this dataset, type `?pbmcsca`

``` r
InstallData("pbmcsca")
data("pbmcsca")
pbmcsca <- NormalizeData(pbmcsca)
pbmcsca <- FindVariableFeatures(pbmcsca)
pbmcsca <- RunFastMNN(object.list = SplitObject(pbmcsca, split.by = "Method"))
pbmcsca <- RunUMAP(pbmcsca, reduction = "mnn", dims = 1:30)
pbmcsca <- FindNeighbors(pbmcsca, reduction = "mnn", dims = 1:30)
pbmcsca <- FindClusters(pbmcsca)
DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/fast_mnn_files/figure-gfm/pbmcsca-1.png)<!-- -->

### Interferon-stimulated and control PBMC

To learn more about this dataset, type `?ifnb`

``` r
InstallData("ifnb")
data("ifnb")
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- RunFastMNN(object.list = SplitObject(ifnb, split.by = "stim"))
ifnb <- RunUMAP(ifnb, reduction = "mnn", dims = 1:30)
ifnb <- FindNeighbors(ifnb, reduction = "mnn", dims = 1:30)
ifnb <- FindClusters(ifnb)
DimPlot(ifnb, group.by = c("stim", "ident", "seurat_annotations"), ncol = 3)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/fast_mnn_files/figure-gfm/ifnb_stim-1.png)<!-- -->

### Eight human pancreatic islet datasets

To learn more about this dataset, type `?panc8`

``` r
InstallData("panc8")
data("panc8")
panc8 <- NormalizeData(panc8)
panc8 <- FindVariableFeatures(panc8)
panc8 <- RunFastMNN(object.list = SplitObject(panc8, split.by = "replicate")[c("celseq", "celseq2", 
    "fluidigmc1", "smartseq2")])
panc8 <- RunUMAP(panc8, reduction = "mnn", dims = 1:30)
panc8 <- FindNeighbors(panc8, reduction = "mnn", dims = 1:30)
panc8 <- FindClusters(panc8)
DimPlot(panc8, group.by = c("replicate", "ident", "celltype"), ncol = 3)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/fast_mnn_files/figure-gfm/pancreas-1.png)<!-- -->
