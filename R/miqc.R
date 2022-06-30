#' @include internal.R
#'
NULL

#' Run miQC on a Seurat object
#'
#' @param object Seurat object
#' @param percent.mt (character) Name of the column in the Seurat metadata that
#'    contains the percent of reads attributed to mitochondrial genes.
#'    Defaults to "percent.mt".
#' @param nFeature_RNA (character) Name of the column in the Seurat metadata that
#'    contains the number of reads per cell. Defaults to "nFeature_RNA".
#' @param posterior.cutoff numeric) The posterior probability of a cell being
#'    part of the compromised distribution, a number between 0 and 1. Any cells
#'    below the appointed cutoff will be marked to keep.
#'    Defaults to 0.75.
#' @param model.type (character) What type of model to generate. A linear
#'    mixture model ("linear") is recommended, but currently b-spline ("spline")
#'    and two-degree polynomial ("polynomial") are also supported
#'    Default = "linear".
#' @param backup.option (character) In case flexmix fails to build a 2 cluster
#'    mixture model, what should RunMiQC do: "percent" (set miQC.keep values
#'    according to backup.percent), "percentile" (set miQC.keep values according
#'    to backup.percentile), "pass" (return original Seurat object), or "halt"
#'    (stop RunMiQC). "percent", "percentile", and "pass" are useful when
#'    processing multiple Seurat objects sequentially.
#' @param backup.percentile (numeric) What percentile to use as cutoff in case
#'    flexmix fails to build a 2 cluster mixture model. Will only be used if
#'    backup.option is "percentile".
#' @param backup.percent (numeric) What percent to use as cutoff in case flexmix
#'    fails to build a 2 cluster mixture model. Will only be used if
#'    backup.option is "percent".
#' @param verbose Boolean. TRUE to show progress messages, FALSE to hide progress messages
#' @details (Copied verbatim from miQC) _Function to fit a two-distribution mixture model on a Seurat object and find those cells probabistically determined to be compromised by the mixture model._
#'
#' @return Returns a Seurat object with probabilities and "keep" decisions stored as "miQC.probability" and "miQC.keep" in the object metadata, respectively.
#' @references Hippen et al. (2021) miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data. bioRxiv doi: 10.1101/2021.03.03.433798
#'
#' @importFrom rlang %||%
#'
#' @export
RunMiQC <- function(
  object,
  percent.mt = "percent.mt",
  nFeature_RNA = "nFeature_RNA",
  posterior.cutoff = 0.75,
  model.type = "linear",
  model.slot = "flexmix_model",
  verbose = TRUE,
  backup.option = "percentile",
  backup.percentile = 0.99,
  backup.percent = 5,
  ...
) {
  SeuratWrappers:::CheckPackage(package = 'flexmix', repository = "CRAN")

  my_data <- Seurat::FetchData(object, vars = c(percent.mt, nFeature_RNA))
  colnames(my_data) <- c("percent.mt", "nFeature_RNA")


  if(!(model.type %in% c("linear", "spline", "polynomial"))){
    stop("model.type must be one of \"linear\", \"spline\", or \"polynomial\"")
  }

  #implementing tryCatch because of internal flexmix error when model fitting
  #fails. see https://github.com/satijalab/seurat-wrappers/issues/108
  my_model <- tryCatch({
    if (model.type == "linear") {
      my_model <- flexmix::flexmix(percent.mt~nFeature_RNA,
                                   data = my_data, k = 2)
    } else if (model.type == "spline") {
      my_model <- flexmix::flexmix(percent.mt~splines::bs(nFeature_RNA),
                                   data = my_data, k = 2)
    } else if (model.type == "polynomial") {
      my_model <- flexmix::flexmix(percent.mt~poly(nFeature_RNA, degree = 2),
                                   data = my_data, k = 2)
    }
  }, error=function(e){
    cat("flexmix fitting error:", conditionMessage(e),"\n")
    my_model <- NULL})

  #
  #set a variable = model_status to denote the status of the model fitting
  #1 = model fit successfully
  #2 = model fit only 1 cluster
  #3 = model fails at flexmix stage
  #

  if(is.null(my_model)){

    model_status <- 3
    model_message <- "flexmix internal failure"
  } else if (ncol(flexmix::parameters(my_model)) == 1){
    model_status <- 2
    model_message <- "flexmix returned only 1 cluster"
  } else if (ncol(flexmix::parameters(my_model)) == 2){
    model_status <- 1
    model_message <- "flexmix model fit successfully"
  } else {
    stop("model_status error, please post issue on GitHub")
  }

  if(model_status %in% c(2,3)){

    if(model_status == 2){
      warning(model_message)
    }
    if(model_status == 3){
      warning(model_message)
    }

    if(backup.option == "halt"){
      stop("Halting.")}
    else if (backup.option == "pass") {
      message("returning object without miQC model or stats")
      return(object)}
    else if (backup.option == "percentile") {
      message("defaulting to backup.percentile for filtering")
      compromised_probability <- 0
      raw_values <- my_data[,percent.mt]
      percentile_cutoff <- quantile(raw_values, probs = backup.percentile)
      cells_to_keep <- ifelse(raw_values <= percentile_cutoff, "keep", "discard")}
    else if (backup.option == "percent"){
      message("defaulting to backup.percent for filtering")
      compromised_probability <- 0
      raw_values <- my_data[,percent.mt]
      cells_to_keep <- ifelse(raw_values <= backup.percent, "keep", "discard")}
    else {
      stop("backup.option must be one of \"percentile\", \"percent\", \"halt\", or \"pass\"")
    }
  } else if (model_status == 1){

    Misc(object, model.slot) <- my_model

    my_model_parameters <- flexmix::parameters(my_model)
    my_model_posterior <- flexmix::posterior(my_model)

    intercept1 <- my_model_parameters[,1][1]
    intercept2 <- my_model_parameters[,2][1]
    if (intercept1 > intercept2) {
      compromised_dist <- 1
    } else {
      compromised_dist <- 2
    }
    compromised_probability <- my_model_posterior[,compromised_dist]
    cells_to_keep <- ifelse(compromised_probability <= posterior.cutoff, "keep", "discard")

  } else {
    stop("model_status error, please post issue on GitHub")
  }

  object <- Seurat::AddMetaData(object = object,
                                metadata = compromised_probability,
                                col.name = "miQC.probability")
  object <- Seurat::AddMetaData(object = object,
                                metadata = cells_to_keep,
                                col.name = "miQC.keep")
  object <- Seurat::LogSeuratCommand(object)

  return(object)
}


#' Run miQC on a Seurat object
#'
#' @param object Seurat object
#' @details _Function to plot the miQC mixture model stored in a Seurat object. `RunMiQC` must be run prior to plotting._
#' @return
#' @references Hippen et al. (2021) miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data. bioRxiv doi: 10.1101/2021.03.03.433798
#'
#' @importFrom rlang %||%
#'
#' @export
PlotMiQC <- function(seurat_object,
                     percent.mt = "percent.mt",
                     nFeature_RNA = "nFeature_RNA",
                     model.slot = "flexmix_model",
                     color.by = "miQC.probability") {

  features_to_fetch <- c(percent.mt, nFeature_RNA, "miQC.probability", "miQC.keep", color.by)
  features_to_fetch <- unique(features_to_fetch)
  my_data <- Seurat::FetchData(seurat_object, vars = features_to_fetch)
  colnames(my_data)[1:2] <- c("percent.mt", "nFeature_RNA")
  my_model <- Misc(seurat_object, model.slot)
  my_model_parameters <- cbind(my_data, flexmix::fitted(my_model))

  #<<< code from plotModel in miQC package >>>
  intercept1 <- flexmix::parameters(my_model, component = 1)[1]
  intercept2 <- flexmix::parameters(my_model, component = 2)[1]
  if (intercept1 > intercept2) {
    compromised_dist <- 1
  } else {
    compromised_dist <- 2
  }

  ggplot2::ggplot(my_data, ggplot2::aes(x = nFeature_RNA, y = percent.mt,
                                        colour = !!ggplot2::sym(color.by))) +
    ggplot2::labs(x = "Unique genes found", y = "Percent reads mitochondrial",
                  color = color.by) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data = my_model_parameters, inherit.aes = FALSE,
                       ggplot2::aes(x = nFeature_RNA, y = Comp.1), lwd = 2) +
    ggplot2::geom_line(data = my_model_parameters, inherit.aes = FALSE,
                       ggplot2::aes(x = nFeature_RNA, y = Comp.2), lwd = 2) +
    ggplot2::ylim(c(0, NA))+
    cowplot::theme_cowplot()
}





