# Load the packages used in the analysis below.
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(fastTopics)

# Set seed to generate results that are reproducible.
set.seed(1)

# Load (and, if necessary, install) the PBMC 3k data.
InstallData("pbmc3k")
data(pbmc3k)

# No pre-processing or pre-selection of genes is needed.

# Fit the multinomial topic model to raw UMI counts. Note that it may
# take several minutes for the model fitting to complete on this data
# set.
pbmc3k <- FitTopicModel(pbmc3k,k = 6)

# This plot shows the cells projected onto the two principal
# components (PCs) of the topic mixture proportions.
Idents(pbmc3k) <- pbmc3k$seurat_annotations
DimPlot(pbmc3k,reduction = "pca_topics",pt.size = 1)

# Compare this against the top two PCs of the transformed counts.
pbmc3k <- FindVariableFeatures(pbmc3k)
pbmc3k <- ScaleData(pbmc3k)
pbmc3k <- RunPCA(pbmc3k)
DimPlot(pbmc3k,reduction = "pca",pt.size = 1)

# Once fitted topic model is extracted, many functions from the
# fastTopics package can be used for analysis and visualization. For
# example, the Structure plot provides an evocative visual summary of
# the estimated mixture proportions for each cell.
fit <- Misc(Reductions(pbmc3k,"multinom_topic_model"))
structure_plot(fit,grouping = Idents(pbmc3k),gap = 25)


