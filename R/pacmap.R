RunPaCMAP <- function(
    object,
    n_components=2,
    n_neighbors=10,
    MN_ratio=0.5,
    FP_ratio=2.0,
    pair_neighbors=None,
    pair_MN=None,
    pair_FP=None,
    distance="euclidean",
    lr=1.0,
    num_iters=(100, 100, 250),
    verbose=False,
    apply_pca=True
){
    if (!py_module_available("pacmap")){
        stop("Please install the pacmap package to use this function.")
    }
    if (!py_module_available("numpy")){
        stop("Please install the numpy package to use this function.")
    }
    if (!py_module_available("annoy")){
        stop("Please install the annoy package to use this function.")
    }
    pacmap_import <- import("pacmap")
    pacmap_import.args <- list(
        n_components=n_components,
        n_neighbors=n_neighbors,
        MN_ratio=MN_ratio,
        FP_ratio=FP_ratio,
        pair_neighbors=pair_neighbors,
        pair_MN=pair_MN,
        pair_FP=pair_FP,
        distance=distance,
        lr=lr,
        num_iters=num_iters,
        verbose=verbose,
        apply_pca=apply_pca
    )
    pacmap <- do.call(pacmap_import$PaCMAP, pacmap_import.args)
    pacmap$fit_transform(as.matrix(x=object))
    return(pacmap$embedding_)

}