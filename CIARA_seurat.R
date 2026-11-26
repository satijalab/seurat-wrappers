#' @include internal.R
#'
NULL

#' get_background_full_seurat
#'
#' @param object A \code{Seurat} object
#' @param threshold threshold in expression for a given gene
#' @param n_cells_low minimum number of cells where a gene is expressed at a
#' level above threshold
#' @param n_cells_high maximum number of cells where a gene is expressed at a
#' level above threshold
#' @return Character vector with all genes expressed at a level higher than
#' \emph{threshold} in a number of cells between \emph{n_cells} and
#' \emph{n_cells_high}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{ https://CRAN.R-project.org/package=CIARA}
#'
#'
#' @export get_background_full_seurat
#'
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom stats fisher.test
#' @importFrom stats median
#' @importFrom Seurat GetAssayData Graphs
get_background_full_seurat <- function(object, threshold = 1, n_cells_low = 3, n_cells_high = 20) {
  CheckPackage(package = "CIARA", repository = "cran")
  if (! (inherits(x= object, what="Seurat"))) {
    stop("Object must be a Seurat object")
  }
  norm_matrix <- as.matrix(GetAssayData(object, slot = "data",assay="RNA"))
  genes_filter <- apply(norm_matrix, 1, function(x) {
    x <- x[x > threshold]
    if (length(x) >= n_cells_low) {
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  })
  
  genes_filter_2 <- apply(norm_matrix, 1, function(x) {
    x <- x[x > threshold]
    if (length(x) > n_cells_high) {
      return(FALSE)
    }
    else{
      return(TRUE)
    }
  })
  genes_important <- row.names(norm_matrix)[genes_filter]
  genes_important_2 <- row.names(norm_matrix)[genes_filter_2]
  genes_important_final <- intersect(genes_important, genes_important_2)
  if (length(genes_important_final) == 0) {
    stop("There are not genes that pass the filtering steps")
  }
  logic_background <- row.names(norm_matrix)
  logic_background[row.names(norm_matrix) %in% genes_important_final] <- TRUE
  logic_background[! (row.names(norm_matrix) %in% genes_important_final)] <- FALSE
  logic_background <- as.logical(logic_background)
  
  object[["RNA"]][['CIARA_background']] <- logic_background
  
  
  return(object)
}










#' CIARA_gene_seurat
#'
#' The gene expression is binarized (1/0) if the value in a given cell is
#' above/below the median. Each of cell with its first K nearest neighbors
#' defined a local region. If there are at least \emph{local_region} enriched
#' in 1 according to \emph{fisher.test}, then the gene is defined as highly
#' localized and a final p value is assigned to it. The final p value is the
#' minimum of the p values from all the enriched local regions. If there are no
#' enriched local regions, then the p value by default is set to 1
#'
#' @param norm_matrix Norm count matrix (n_genes X n_cells)
#' @param knn_matrix K-nearest neighbors matrix (n_cells X n_cells).
#' @param gene_expression numeric vector with the gene expression (length equal
#' to n_cells). The gene expression is binarized (equal to 0/1 in the cells
#' where the value is below/above the median)
#' @param p_value p value returned by the function \emph{fisher.test} with
#' parameter alternative = "g"
#' @param odds_ratio odds_ratio returned by the function \emph{fisher.test}
#' with parameter alternative = "g"
#' @param local_region Integer. Minimum number of local regions (cell with its
#' knn neighbours) where the binarized gene expression is enriched in 1.
#' @param approximation Logical.For a given gene, the fisher test is run in the
#' local regions of only the cells where the binarized gene expression is 1.
#' @inheritParams get_background_full_seurat
#' @return List with one element corresponding to the p value of the gene.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{ https://CRAN.R-project.org/package=CIARA}
#' @export CIARA_gene_seurat
CIARA_gene_seurat <- function(norm_matrix, knn_matrix, gene_expression, p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = FALSE) {
  
  
  
  
  median_genes <- median(gene_expression)
  binary_expression <- rep(0, length(gene_expression))
  binary_expression[gene_expression > median_genes] <- 1
  binary_expression[gene_expression <= median_genes] <- 0
  
  
  if (approximation == FALSE) {
    sub_feature <- colnames(norm_matrix)
    message("approximation == FALSE")
  }
  else {
    sub_feature <- colnames(norm_matrix)[which(binary_expression == 1)]
    message("approximation == TRUE")
  }
  if (!all(colnames(norm_matrix)%in%row.names(knn_matrix))) {
    stop ("row.names in knn_matrix are not equal to colnames in norm_matrix")
    
  }
  
  knn_matrix_big <- knn_matrix[sub_feature, ]
  
  
  
  perform_fisher <- function(j) {
    
    nn_cell <- colnames(norm_matrix)[which(knn_matrix_big[sub_feature[j], ] > 0)]
    nn_gene_expression <- binary_expression[colnames(norm_matrix)%in%nn_cell]
    if(sum(nn_gene_expression) != 0) {
      input_fisher <- matrix(c(sum(nn_gene_expression == 1), sum(binary_expression == 1)-sum(nn_gene_expression == 1), sum(nn_gene_expression == 0), sum(binary_expression == 0)-sum(nn_gene_expression == 0)), nrow = 2, ncol = 2, byrow = T)
      test_output <- fisher.test(input_fisher, alternative = "g")
      return(test_output$p.value)
    }
    else{
      return(1)
    }
  }
  
  p_value_gene <- unlist(lapply(seq_len(length(sub_feature)), perform_fisher))
  
  if(sum(p_value_gene < p_value) >= local_region) {
    p_value_final <- min(p_value_gene)
  }
  else {
    p_value_final <- 1
  }
  if (all(p_value_final == 1)){
    warning(paste0("There are not genes enriched in ", local_region, " or more local regions"))
  }
  
  return((p_value_final))
}




#' CIARA_seurat
#'
#' It selects highly localized genes as specified in \emph{CIARA_gene_seurat},
#' starting from genes in \emph{background}
#'
#' @param background Vector of genes for which the function \emph{CIARA_gene_seurat}
#' is run.
#' @param cores_number Integer.Number of cores to use.
#' @inheritParams CIARA_gene_seurat
#' @inheritParams get_background_full_seurat
#' @return Dataframe with n_rows equal to the length of
#' \emph{background} . Each row is the output from \emph{CIARA_gene_seurat}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{ https://CRAN.R-project.org/package=CIARA}
#'
#'
#' @export CIARA_seurat
CIARA_seurat <- function(object, background, cores_number = 1, p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = FALSE) {
  CheckPackage(package = "CIARA", repository = "cran")
  if (! (inherits(x= object, what="Seurat"))) {
    stop("Object must be a Seurat object")
  }
  norm_matrix <- as.matrix(GetAssayData(object, slot = "data",assay="RNA"))
  
  if (! ("RNA_nn" %in% Graphs(object))){
    stop ("RNA_nn slot not present in seurat object. Run first the function FindNeighbors")
  }
  knn_matrix <- as.matrix(object[["RNA_nn"]])
  if (! ("CIARA_background" %in% colnames(object[["RNA"]][[]]))) {
    stop("Run first the function get_background_full_seurat before CIARA")
  }
  background <- row.names(object[["RNA"]][['CIARA_background']])[object[["RNA"]][['CIARA_background']][,1]]
  run_loop_genes = function(i) {
    if (!all(background %in% row.names(norm_matrix))) {
      stop("Some background genes are not present in norm matrix")
    }
    gene_expression <- as.vector(norm_matrix[background[i],])
    message(paste0("Running CIARA on gene:", background[i]))
    return(CIARA_gene_seurat(norm_matrix, knn_matrix, gene_expression, p_value, odds_ratio, local_region, approximation))
  }
  result_final <- do.call(rbind, mclapply(seq_len(length(background)),run_loop_genes , mc.cores = cores_number))
  names(result_final)=background
  object[["RNA"]][['CIARA_p_value']] <- result_final
  
  return(object)
}



















