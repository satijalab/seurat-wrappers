# Load the packages used in the analysis below.
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(fastTopics)

# Set seed to generate results that are reproducible.

# Load the UMI count data.
InstallData("pbmc3k")
data(pbmc3k)

# No pre-processing is needed.

# Fit topic model to UMI counts. Note that it may take several minutes
# for the model fitting to complete.
pbmc3k <- FitTopicModel(pbmc3k,k = 6,numiter.main = 10,numiter.refine = 10)

# Plot the top two PCs of the mixture proportions.
Idents(pbmc3k) <- pbmc3k$seurat_annotations
DimPlot(pbmc3k)

# Create a Structure plot.
fit <- Misc(Reductions(pbmc3k,"multinom_topic_model"))
p <- structure_plot(fit,grouping = Idents(pbmc3k),gap = 30)

# Fit non-negative matrix factorization to UMI counts.
# TO DO.

