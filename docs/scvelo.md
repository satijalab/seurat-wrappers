Estimating RNA Velocity using Seurat and scVelo
================
Compiled: May 26, 2020

This vignette demonstrates analysing RNA Velocity quantifications stored
in a Seurat object.
<!-- Parameters are based off of the [RNA Velocity tutorial](http://pklab.med.harvard.edu/velocyto/notebooks/R/SCG71.nb.html).-->
If you use scVelo in your work, please
cite:

<!-- > *RNA velocity of single cells* -->

<!-- > -->

<!-- > Gioele La Manno, Ruslan Soldatov, Amit Zeisel, Emelie Braun, Hannah Hochgerner, Viktor Petukhov, Katja Lidschreiber, Maria E. Kastriti, Peter Lönnerberg, Alessandro Furlan, Jean Fan, Lars E. Borm, Zehua Liu, David van Bruggen, Jimin Guo, Xiaoling He, Roger Barker, Erik Sundström, Gonçalo Castelo-Branco, Patrick Cramer, Igor Adameyko, Sten Linnarsson & Peter V. Kharchenko -->

<!-- > -->

<!-- > doi: [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6) -->

<!-- > -->

<!-- > Website: https://velocyto.org -->

Prerequisites to install:

  - [Seurat](https://satijalab.org/seurat/install)
  - [scVelo](https://scvelo.readthedocs.io/)
  - [SeuratDisk](https://mojaveazure.github.io/seurat-disk)
  - [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)

<!-- end list -->

``` r
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
```

``` r
# If you don't have velocyto's example mouse bone marrow dataset, download with the CURL command
# curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile
# = '~/Downloads/SCG71.loom')
ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
bm <- as.Seurat(x = ldat)
bm[["RNA"]] <- bm[["spliced"]]
DefaultAssay(bm) <- "RNA"
bm <- NormalizeData(bm)
bm <- FindVariableFeatures(bm)
bm <- ScaleData(bm, features = rownames(bm))
bm <- RunPCA(bm, features = VariableFeatures(bm))
bm <- FindNeighbors(bm)
bm <- FindClusters(bm)
bm <- RunUMAP(bm, dims = 1:20)
SaveH5Seurat(bm, filename = "mouseBM.h5Seurat")
Convert("mouseBM.h5Seurat", dest = "h5ad")
```

In Python

``` python
import scvelo
adata = scvelo.read("mouseBM.h5ad")
adata
```

    ## AnnData object with n_obs × n_vars = 6667 × 24421
    ##     obs: 'orig.ident', 'nCount_spliced', 'nFeature_spliced', 'nCount_unspliced', 'nFeature_unspliced', 'nCount_ambiguous', 'nFeature_ambiguous', 'nCount_RNA', 'nFeature_RNA', 'RNA_snn_res.0.8', 'seurat_clusters'
    ##     var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable', 'ambiguous_features', 'spliced_features', 'unspliced_features'
    ##     uns: 'neighbors'
    ##     obsm: 'X_pca', 'X_umap'
    ##     varm: 'PCs'
    ##     layers: 'ambiguous', 'spliced', 'unspliced'
    ##     obsp: 'distances'

``` python
scvelo.pp.filter_and_normalize(adata)
```

    ## Normalized count data: spliced, unspliced.
    ## WARNING: Did not modify X as it looks preprocessed already.

``` python
scvelo.pp.moments(adata)
```

    ## computing neighbors
    ##     finished (0:00:04) --> added 
    ##     'distances' and 'connectivities', weighted adjacency matrices (adata.obsp)
    ## computing moments based on connectivities
    ##     finished (0:00:04) --> added 
    ##     'Ms' and 'Mu', moments of spliced/unspliced abundances (adata.layers)

``` python
scvelo.tl.velocity(adata)
```

    ## computing velocities
    ##     finished (0:00:11) --> added 
    ##     'velocity', velocity vectors for each individual cell (adata.layers)

``` python
scvelo.tl.velocity_graph(adata)
```

    ## computing velocity graph
    ## ... 13%... 27%... 41%... 54%... 68%... 81%... 95%... 100%    finished (0:00:22) --> added 
    ##     'velocity_graph', sparse matrix with cosine correlations (adata.uns)

``` python
scvelo.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters")
```

    ## computing velocity embedding
    ##     finished (0:00:01) --> added
    ##     'velocity_umap', embedded velocity vectors (adata.obsm)

<img src="scvelo_files/figure-gfm/scvelo-1.png" width="1536" />

``` python
scvelo.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters")
```

<img src="scvelo_files/figure-gfm/scvelo-2.png" width="1536" />
