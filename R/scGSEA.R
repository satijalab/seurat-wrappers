#' @include internal.R
#'
NULL

#' Run Single cell Gene Set Enrichement Analysis on GF-ICF on a Seurat object
#'
#' This function run runScGSEA function from GF-ICF package on Seurat object. 
#' It computes GSEA for each cells across a set of input pathways by using NMF.
#'
#' @param object Seurat object
#' @param assay Assay to use, defaults to the default assay
#' @param filterGenes Rank of NMF to use for the enrichment
#' @param nmf.dim Rank of NMF to use for the enrichment
#' @param geneID The type of gene names in rownames of \code{object}. It can be either 'ensamble' or 'symbol'. Default: 'ensamble'
#' @param species The type of species in \code{object}. It can be either 'mouse' or 'human'. Default: 'human'
#' @param category MSigDB collection abbreviation, such as H or C1
#' @param subcategory MSigDB sub-collection abbreviation, such as CGP or BP
#' @param rescale If different by none, pathway's activity scores are resealed as Z-score. Possible values are none, byGS or byCell. Default is none
#'
#' @return returns a Seurat object with pathways x cell matrix 
#' @import RcppML
#' @export

RunScGSEA <- function(
    object,
    assay = NULL,
    filterGenes = TRUE,
    nmf.dim = 100,
    geneID = c("ensamble","symbol"),
    species = c("mouse", "human"),
    category = "H",
    subcategory = NULL,
    verbose = F) {
      
    SeuratWrappers:::CheckPackage(package = 'gambalab/gficf', repository = 'github')
    assay <- assay %||% DefaultAssay(object = object)
    
    # Get raw count matrix from Seurat object
    M <- GetAssayData(object = object, slot = "counts")
    # Data normalization and gene filtering
    M <- gficf::gficf(M = M, filterGenes = filterGenes)
    # Create NMF-subspace 
    M <- gficf::runNMF(data = M, dim = nmf.dim)
    # Create t-UMAP space
    M <- gficf::runReduction(data = M, reduction = "umap", verbose = verbose)
    # Compute GSEA for each cells across a set of input pathways by using NMF
    M <- gficf::runScGSEA(data = M,
                          geneID = geneID,
                          species = species,
                          category = category,
                          subcategory = subcategory,
                          nmf.k = nmf.dim,
                          fdr.th = .1,
                          rescale = 'none',
                          verbose = verbose)
    
    # Add cell names to enrichment results
    raw_enrich <- t(M$scgsea$x)
    colnames(raw_enrich) <- colnames(object)
    # Normalize enrichment results by computing pathway's activity z-scores
    norm_enrich <- Matrix::Matrix((raw_enrich - Matrix::rowMeans(raw_enrich)) / apply(raw_enrich, 1, sd), sparse=T)
    # Create a new Seurat object containing the raw pathway x cell matrix in counts slot and preserving meta data (adding _byGenes tag only if clusters were already computed)
    if('seurat_clusters' %in% colnames(object@meta.data)) {
      colnames(object@meta.data) <- gsub('seurat_clusters', 'seurat_clusters_byGenes', colnames(object@meta.data))
    }
    path_obj <- CreateSeuratObject(counts = raw_enrich, meta.data = object@meta.data)
    # And the z-score-normalized one in data and scale.data slots
    path_obj <- SetAssayData(object = path_obj, slot = 'data', new.data = norm_enrich)
    path_obj <- SetAssayData(object = path_obj, slot = 'scale.data', new.data = as.matrix(norm_enrich))
    # Store metadata of pathways retained after significance filtering
    feat.meta <- M$scgsea$stat[M$scgsea$stat$pathway%in%colnames(M$scgsea$x), ] 
    feat.meta <- data.frame(feat.meta, row.names = 1)
    feat.meta$genes <- do.call(c, lapply(M$scgsea$pathways[rownames(feat.meta)], function(x) paste(x, collapse = ',')))
    path_obj[[assay]]@meta.features <- feat.meta
    
    # Store dimensionality reduction results computed on genes x cells matrix
    path_obj@reductions <- object@reductions
    names(path_obj@reductions) <- paste0(names(path_obj@reductions), '_byGenes')
    
    return(path_obj)
      
}


