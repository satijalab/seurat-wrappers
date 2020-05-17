Import alevin counts & generate Seurat object
================
Compiled: May 16, 2020

This vignette demonstrates the import of alevin quantified counts into
Seurat. Commands and parameters are based off of the [alevin
tutorial](https://combine-lab.github.io/alevin-tutorial/2018/running-alevin/).
If you use alevin in your work, please cite:

> *Alevin efficiently estimates accurate gene abundances from dscRNA-seq
> data*
> 
> Avi Srivastava, Laraib Malik, Tom Smith, Ian Sudbery & Rob Patro
> 
> Genome Biology, 2019.
> 
> doi:
> [10.1186/s13059-019-1670-y](https://doi.org/10.1186/s13059-019-1670-y)
> 
> GitHub: <https://github.com/COMBINE-lab/salmon>

Prerequisites to install:

  - [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
  - [AlevinRtools](https://github.com/k3yavi/alevin-Rtools)

<!-- end list -->

``` r
library(SeuratWrappers)
library(AlevinRtools)
```

## 

### Import alevin quantified counts

``` r
pbmc <- ReadAlevin("~/alevin_out/alevin/quants_mat.gz")
```
