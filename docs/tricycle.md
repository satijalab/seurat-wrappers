Running estimate\_cycle\_position from tricycle on Seurat Objects
================
Compiled: July 09, 2021

  - [Introduction](#introduction)
  - [Loading examle data and making Seurat
    object](#loading-examle-data-and-making-seurat-object)
  - [Inferring the cell cycle
    position](#inferring-the-cell-cycle-position)
  - [Visualizing the results](#visualizing-the-results)
  - [Assessing performance](#assessing-performance)
  - [Plot out the kernel density](#plot-out-the-kernel-density)
  - [Resoures for tricycle](#resoures-for-tricycle)

This vignette demonstrates the use of the estimate\_cycle\_position from
the tricycle package on Seurat objects.

> *Universal prediction of cell cycle position using transfer learning*
> 
> Shijie C. Zheng, Genevieve Stein-O’Brien, Jonathan J. Augustin, Jared
> Slosberg, Giovanni A. Carosso, Briana Winer, Gloria Shin, Hans T.
> Bjornsson, Loyal A. Goff, Kasper D. Hansen
> 
> bioRxiv, 2021.
> 
> doi:
> [10.1101/2021.04.06.438463](https://doi.org/10.1101/2021.04.06.438463)
> 
> Bioconductor:
> <https://www.bioconductor.org/packages/release/bioc/html/tricycle.html>

Prerequisites to install:

  - [Seurat](https://satijalab.org/seurat/install)
  - [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
  - [tricycle](https://www.bioconductor.org/packages/release/bioc/html/tricycle.html)

<!-- end list -->

``` r
library(Seurat)
library(SeuratWrappers)
library(tricycle)
```

## Introduction

The Biocondutor package
[tricycle](https://www.bioconductor.org/packages/release/bioc/html/tricycle.html)
infers cell cycle position for a single-cell RNA-seq dataset. Here, we
show the implementation of **main** function of tricycle,
estimate\_cycle\_position, on the Seurat objects. More information can
be found at
[tricycle](https://www.bioconductor.org/packages/release/bioc/html/tricycle.html).

## Loading examle data and making Seurat object

``` r
data(neurosphere_example, package = "tricycle")
neurosphere_example <- as.Seurat(neurosphere_example)
neurosphere_example
```

    ## An object of class Seurat 
    ## 1500 features across 400 samples within 1 assay 
    ## Active assay: RNA (1500 features, 0 variable features)

Note that after converting the SingleCellExperiment object to Seurat
object, the original “logcounts” assay is saved as a slot with name
“data” in Seurat default Assay.

## Inferring the cell cycle position

The `Runtricycle()` function in the SeuratWrappers package first project
the data into the cell cycle embeddings using the internal reference in
tricycle package, and then estimate the cell cycle position. The
estimated cell cycle position is bound between 0 and 2pi. Note that we
strive to get high resolution of cell cycle state, and we think the
continuous position is more appropriate when describing the cell cycle.
However, to help users understand the position variable, we also note
that users can approximately relate 0.5pi to be the start of S stage, pi
to be the start of G2M stage, 1.5pi to be the middle of M stage, and
1.75pi-0.25pi to be G1/G0 stage.

``` r
neurosphere_example <- Runtricycle(object = neurosphere_example, slot = "data", reduction.name = "tricycleEmbedding", 
    reduction.key = "tricycleEmbedding_", gname = NULL, gname.type = "ENSEMBL", species = "mouse")
```

## Visualizing the results

We could extract the cell cycle embedding and make a scatter plot of the
embeddings colored by the position inference. And we also extract the
expression level of gene Top2a for accessing the performance, described
below.

``` r
plot.df <- FetchData(object = neurosphere_example, vars = c("tricycleEmbedding_1", "tricycleEmbedding_2", 
    "tricyclePosition", "ENSMUSG00000020914"))
names(plot.df)[4] <- "Top2a"
```

Let us plot out the cell cycle embedding. You could also plot other
embeddings, such as T\_SNE or UMAP with points colored by the cell cycle
position.

``` r
library(ggplot2)
library(cowplot)
p <- tricycle:::.plot_emb_circle_scale(emb.m = plot.df[, 1:2], color.value = plot.df$tricyclePosition, 
    color_by = "tricyclePosition", point.size = 3.5, point.alpha = 0.9)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/tricycle_files/figure-gfm/plotemb-1.png)<!-- -->

## Assessing performance

We have two ways of (quickly) assessing whether triCycle works. They are

1.  Look at the projection of the data into the cell cycle embedding.
2.  Look at the expression of key genes as a function of cell cycle
    position.

Plotting the projection of the data into the cell cycle embedding is
shown above. Our observation is that deeper sequenced data will have a
more clearly ellipsoid pattern with an empty interior. As sequencing
depth decreases, the radius of the ellipsoid decreases until the empty
interior disappears. So the absence of an interior does not mean the
method does not work.

It is more important to inspect a couple of genes as a function of cell
cycle position. We tend to use Top2a which is highly expressed and
therefore “plottable” in every dataset. Other candidates are for example
Smc2. To plot this data, we provide a convenient function
`fit_periodic_loess()` to fit a loess line between the cyclic variable
\(\theta\) and other response variables. This fitting is done by making
`theta.v` 3 periods `(c(theta.v - 2 * pi, theta.v, theta.v + 2 * pi))`
and repeating `y` 3 times. Only the fitted values corresponding to
original `theta.v` will be returned. In this example, we show how well
the expression of the cell cycle marker gene *Top2a* change along
\(\theta\).

``` r
fit.l <- fit_periodic_loess(neurosphere_example$tricyclePosition, plot.df$Top2a, plot = TRUE, x_lab = "Cell cycle position θ", 
    y_lab = "log2(Top2a)", fig.title = paste0("Expression of Top2a along θ (n=", ncol(neurosphere_example), 
        ")"))
names(fit.l)
```

    ## [1] "fitted"   "residual" "pred.df"  "loess.o"  "rsquared" "fig"

``` r
fit.l$fig + theme_bw(base_size = 14)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/tricycle_files/figure-gfm/loess-1.png)<!-- -->

For Top2a we expect peak expression around \(\pi\).

## Plot out the kernel density

Another useful function is *plot\_ccposition\_den*, which computes
kernel density of \(\theta\) conditioned on a phenotype using von Mises
distribution. The ouput figures are provided in two flavors, polar
coordinates and Cartesian coordinates. This could be useful when
comparing different cell types, treatments, or just stages. (Because we
use a very small dataset here as example, we set the bandwith, i.e. the
concentration parameter of the von Mises distribution as 10 to get a
smooth line.)

``` r
plot_ccposition_den(neurosphere_example$tricyclePosition, neurosphere_example$sample, "sample", 
    bw = 10, fig.title = "Kernel density of θ") + theme_bw(base_size = 14)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/tricycle_files/figure-gfm/density-1.png)<!-- -->

``` r
plot_ccposition_den(neurosphere_example$tricyclePosition, neurosphere_example$sample, "sample", 
    type = "circular", bw = 10, fig.title = "Kernel density of θ") + theme_bw(base_size = 14)
```

![](/__w/seurat-wrappers/seurat-wrappers/test-build/tricycle_files/figure-gfm/density-2.png)<!-- -->

## Resoures for tricycle

More information about constructing your own reference, other usages and
running tricycle outside of the Seurat environment can be found at
[tricycle](https://www.bioconductor.org/packages/release/bioc/html/tricycle.html).
