Analyzing Seurat data using fastTopics
================
Compiled: January 15, 2021

*Add introductory text here.*

*Point to the fastTopics website and vignettes.*

*Add citation here.*

We begin by loading the packages used to perform the analysis.

``` r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(fastTopics)
```

We set the seed so that the results can be reproduced.

``` r
set.seed(1)
```

Load—and, if necessary, install—the PBMC 3k data set containing
transcription profiles for 2,700 cells.

``` r
InstallData("pbmc3k")
data(pbmc3k)
dim(GetAssayData(pbmc3k))
# [1] 13714  2700
```

Fit the multinomial topic model to the raw UMI counts—*no pre-processing
or pre-selection of genes is needed.* Note that it may take several
minutes to complete this model fitting step.

``` r
pbmc3k <- FitTopicModel(pbmc3k,k = 6)
```

To fit a topic model, we must specify \(K\), the number of topics. Here,
we have chosen \(K = 6\) topics. In most settings, a good choice of
\(K\) will not be known in advance, so you will you want to explore the
results from topic models at different settings of \(K\).

This plot shows the cells projected onto the top two principal
components (PCs) obtained from the topic mixture proportions.

``` r
Idents(pbmc3k) <- pbmc3k$seurat_annotations
DimPlot(pbmc3k,reduction = "pca_topics",pt.size = 1)
```

<img src="fasttopics_files/figure-gfm/pca-1-1.png" style="display: block; margin: auto;" />

Compare this against PCs 1 and 2 of the transformed counts:

``` r
pbmc3k <- FindVariableFeatures(pbmc3k)
pbmc3k <- ScaleData(pbmc3k)
# Centering and scaling data matrix
pbmc3k <- RunPCA(pbmc3k)
# PC_ 1 
# Positive:  FTL, FTH1, COTL1, CST3, OAZ1, ACTB, LGALS1, S100A4, AIF1, FCER1G 
#      TMSB4X, S100A6, TYROBP, LST1, PSAP, TYMP, SAT1, S100A11, CTSS, SPI1 
#      SERPINA1, LYZ, TMSB10, IFITM3, HLA-DRB1, FCN1, CFD, HLA-DPA1, VIM, GSTP1 
# Negative:  MALAT1, IL32, LTB, CCL5, CTSW, CD247, CD2, NKG7, LINC00926, ACAP1 
#      TCL1A, CST7, GZMA, FGFBP2, NCR3, BEX2, GZMK, HOPX, GNLY, MAL 
#      SAMD3, STK17A, MYC, XCL1, SPON2, NELL2, LDLRAP1, ZAP70, XCL2, PRF1 
# PC_ 2 
# Positive:  FTL, TYROBP, S100A8, S100A9, FCN1, AIF1, FTH1, LYZ, CTSS, LST1 
#      CFD, S100A6, TYMP, SAT1, SERPINA1, LGALS2, S100A11, PSAP, S100A4, CST3 
#      FCER1G, IFITM3, CFP, APOBEC3A, LGALS1, SPI1, LGALS3, TIMP1, NPC2, IFI30 
# Negative:  ACTG1, STMN1, TUBA1B, TYMS, ZWINT, GZMA, KIAA0101, TK1, RRM2, HMGB2 
#      H2AFZ, DUT, HMGA1, FEN1, MYBL2, BIRC5, GINS2, GAPDH, SRSF3, NKG7 
#      HSP90AA1, PTTG1, IL32, FABP5, ASF1B, MKI67, RANBP1, KIFC1, CENPM, ACTB 
# PC_ 3 
# Positive:  CD74, HLA-DRA, HLA-DPB1, HLA-DQB1, HLA-DQA1, HLA-DRB1, HLA-DPA1, HLA-DQA2, HLA-DRB5, LYZ 
#      HLA-DMB, CD1C, CST3, VIM, LGALS2, HLA-DMA, MS4A1, CD79A, CLEC10A, LINC00926 
#      CD79B, GSTP1, MS4A6A, TCL1A, S100A9, FCER1A, PLD4, FCN1, S100A8, IFI6 
# Negative:  PPBP, GNG11, SPARC, PF4, AP001189.4, SDPR, ITGA2B, CLU, CD9, TREML1 
#      GP9, LY6G6F, CMTM5, NRGN, TUBB1, RP11-367G6.3, C6orf25, GP1BA, RGS18, F13A1 
#      CA2, SCGB1C1, RUFY1, CLDN5, SEPT5, CLEC1B, ITGB3, HIST1H2AC, TUBA4A, CTTN 
# PC_ 4 
# Positive:  CD74, HLA-DQB1, HLA-DQA1, HLA-DRA, HLA-DQA2, HLA-DPB1, CD79A, HLA-DRB1, MS4A1, HLA-DPA1 
#      HLA-DRB5, CD79B, HLA-DMB, CD1C, LINC00926, GPX1, TCL1A, HLA-DMA, KIAA0125, PPBP 
#      CLEC10A, P2RX5, PPP1R14A, PF4, PLD4, SPARC, FCER1A, GNG11, ITM2C, ITGA2B 
# Negative:  NKG7, GZMA, GNLY, PRF1, FGFBP2, SPON2, CTSW, FCGR3A, CST7, GZMB 
#      AKR1C3, CLIC3, GZMH, CCL5, CD247, XCL2, HOPX, IFITM2, CCL4, XCL1 
#      FCER1G, TTC38, S100A4, CTSC, IL32, LST1, MS4A7, SERPINA1, AIF1, SAMD3 
# PC_ 5 
# Positive:  ZWINT, KIAA0101, RRM2, HMGB2, AQP3, FEN1, TYMS, H2AFZ, GINS2, KIFC1 
#      TK1, MYBL2, HN1, BIRC5, LTB, S100A8, FTL, NUSAP1, CDC45, CKS1B 
#      DUT, CENPM, MKI67, STMN1, CDCA5, ASF1B, S100A9, GAPDH, CFD, AIF1 
# Negative:  NKG7, CD74, HLA-DQA1, GNLY, SPON2, HLA-DPB1, FGFBP2, HLA-DQA2, HLA-DRA, GZMB 
#      HLA-DQB1, PRF1, HLA-DPA1, HLA-DRB1, GZMA, CTSW, HLA-DRB5, CLIC3, AKR1C3, CST7 
#      CD1C, XCL2, CCL4, XCL1, HOPX, CCL5, CLEC10A, HLA-DMB, TTC38, FCER1A
DimPlot(pbmc3k,reduction = "pca",pt.size = 1)
```

<img src="fasttopics_files/figure-gfm/pca-2-1.png" style="display: block; margin: auto;" />

Once fitted topic model is extracted, many functions from the fastTopics
package can be used for analysis and visualization. For example, the
Structure plot provides an evocative visual summary of the estimated
mixture proportions for each cell.

``` r
fit <- Misc(Reductions(pbmc3k,"multinom_topic_model"))
structure_plot(fit,grouping = Idents(pbmc3k),gap = 25)
# Perplexity automatically changed to 83 because original setting of 100 was too large for the number of samples (255)
# Perplexity automatically changed to 65 because original setting of 100 was too large for the number of samples (201)
# Perplexity automatically changed to 38 because original setting of 100 was too large for the number of samples (120)
# Perplexity automatically changed to 37 because original setting of 100 was too large for the number of samples (115)
# Perplexity automatically changed to 6 because original setting of 100 was too large for the number of samples (24)
# Perplexity automatically changed to 2 because original setting of 100 was too large for the number of samples (10)
```

<img src="fasttopics_files/figure-gfm/structure-plot-1.png" style="display: block; margin: auto;" />

    # Read the 516 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 100.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.15 seconds (sparsity = 0.725084)!
    # Learning embedding...
    # Iteration 50: error is 46.544798 (50 iterations in 0.07 seconds)
    # Iteration 100: error is 44.191974 (50 iterations in 0.06 seconds)
    # Iteration 150: error is 44.191098 (50 iterations in 0.06 seconds)
    # Iteration 200: error is 44.191097 (50 iterations in 0.06 seconds)
    # Iteration 250: error is 44.191097 (50 iterations in 0.05 seconds)
    # Iteration 300: error is 0.267504 (50 iterations in 0.06 seconds)
    # Iteration 350: error is 0.265744 (50 iterations in 0.06 seconds)
    # Iteration 400: error is 0.265741 (50 iterations in 0.06 seconds)
    # Iteration 450: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 500: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 550: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 600: error is 0.265742 (50 iterations in 0.05 seconds)
    # Iteration 650: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 700: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 750: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 800: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 850: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 900: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 950: error is 0.265742 (50 iterations in 0.06 seconds)
    # Iteration 1000: error is 0.265742 (50 iterations in 0.06 seconds)
    # Fitting performed in 1.14 seconds.
    # Read the 358 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 100.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.10 seconds (sparsity = 0.948269)!
    # Learning embedding...
    # Iteration 50: error is 42.965050 (50 iterations in 0.04 seconds)
    # Iteration 100: error is 41.889224 (50 iterations in 0.04 seconds)
    # Iteration 150: error is 41.888879 (50 iterations in 0.04 seconds)
    # Iteration 200: error is 41.888918 (50 iterations in 0.04 seconds)
    # Iteration 250: error is 41.888924 (50 iterations in 0.04 seconds)
    # Iteration 300: error is 0.154191 (50 iterations in 0.04 seconds)
    # Iteration 350: error is 0.154156 (50 iterations in 0.04 seconds)
    # Iteration 400: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 450: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 500: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 550: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 600: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 650: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 700: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 750: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 800: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 850: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 900: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 950: error is 0.154155 (50 iterations in 0.04 seconds)
    # Iteration 1000: error is 0.154155 (50 iterations in 0.04 seconds)
    # Fitting performed in 0.73 seconds.
    # Read the 356 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 100.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.10 seconds (sparsity = 0.948886)!
    # Learning embedding...
    # Iteration 50: error is 44.127036 (50 iterations in 0.04 seconds)
    # Iteration 100: error is 42.960180 (50 iterations in 0.04 seconds)
    # Iteration 150: error is 42.946318 (50 iterations in 0.03 seconds)
    # Iteration 200: error is 42.946371 (50 iterations in 0.04 seconds)
    # Iteration 250: error is 42.946352 (50 iterations in 0.04 seconds)
    # Iteration 300: error is 0.264561 (50 iterations in 0.04 seconds)
    # Iteration 350: error is 0.264515 (50 iterations in 0.03 seconds)
    # Iteration 400: error is 0.264520 (50 iterations in 0.04 seconds)
    # Iteration 450: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 500: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 550: error is 0.264520 (50 iterations in 0.04 seconds)
    # Iteration 600: error is 0.264520 (50 iterations in 0.04 seconds)
    # Iteration 650: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 700: error is 0.264520 (50 iterations in 0.04 seconds)
    # Iteration 750: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 800: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 850: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 900: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 950: error is 0.264520 (50 iterations in 0.03 seconds)
    # Iteration 1000: error is 0.264520 (50 iterations in 0.03 seconds)
    # Fitting performed in 0.71 seconds.
    # Read the 255 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 83.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.05 seconds (sparsity = 0.995309)!
    # Learning embedding...
    # Iteration 50: error is 41.894862 (50 iterations in 0.03 seconds)
    # Iteration 100: error is 41.348727 (50 iterations in 0.02 seconds)
    # Iteration 150: error is 41.316669 (50 iterations in 0.02 seconds)
    # Iteration 200: error is 41.331596 (50 iterations in 0.02 seconds)
    # Iteration 250: error is 41.337796 (50 iterations in 0.02 seconds)
    # Iteration 300: error is 0.153763 (50 iterations in 0.02 seconds)
    # Iteration 350: error is 0.152951 (50 iterations in 0.02 seconds)
    # Iteration 400: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 450: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 500: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 550: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 600: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 650: error is 0.152947 (50 iterations in 0.02 seconds)
    # Iteration 700: error is 0.152948 (50 iterations in 0.02 seconds)
    # Iteration 750: error is 0.152947 (50 iterations in 0.02 seconds)
    # Iteration 800: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 850: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 900: error is 0.152947 (50 iterations in 0.02 seconds)
    # Iteration 950: error is 0.152945 (50 iterations in 0.02 seconds)
    # Iteration 1000: error is 0.152947 (50 iterations in 0.02 seconds)
    # Fitting performed in 0.46 seconds.
    # Read the 201 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 65.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.03 seconds (sparsity = 0.993688)!
    # Learning embedding...
    # Iteration 50: error is 43.455428 (50 iterations in 0.02 seconds)
    # Iteration 100: error is 42.659459 (50 iterations in 0.02 seconds)
    # Iteration 150: error is 42.729750 (50 iterations in 0.02 seconds)
    # Iteration 200: error is 43.218125 (50 iterations in 0.02 seconds)
    # Iteration 250: error is 43.736809 (50 iterations in 0.02 seconds)
    # Iteration 300: error is 0.664267 (50 iterations in 0.02 seconds)
    # Iteration 350: error is 0.649046 (50 iterations in 0.02 seconds)
    # Iteration 400: error is 0.581421 (50 iterations in 0.02 seconds)
    # Iteration 450: error is 0.581122 (50 iterations in 0.02 seconds)
    # Iteration 500: error is 0.581124 (50 iterations in 0.02 seconds)
    # Iteration 550: error is 0.581123 (50 iterations in 0.02 seconds)
    # Iteration 600: error is 0.581120 (50 iterations in 0.02 seconds)
    # Iteration 650: error is 0.581120 (50 iterations in 0.02 seconds)
    # Iteration 700: error is 0.581120 (50 iterations in 0.02 seconds)
    # Iteration 750: error is 0.581124 (50 iterations in 0.02 seconds)
    # Iteration 800: error is 0.581124 (50 iterations in 0.02 seconds)
    # Iteration 850: error is 0.581123 (50 iterations in 0.02 seconds)
    # Iteration 900: error is 0.581124 (50 iterations in 0.02 seconds)
    # Iteration 950: error is 0.581120 (50 iterations in 0.02 seconds)
    # Iteration 1000: error is 0.581120 (50 iterations in 0.02 seconds)
    # Fitting performed in 0.35 seconds.
    # Read the 120 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 38.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.01 seconds (sparsity = 0.988750)!
    # Learning embedding...
    # Iteration 50: error is 45.952359 (50 iterations in 0.01 seconds)
    # Iteration 100: error is 48.966439 (50 iterations in 0.01 seconds)
    # Iteration 150: error is 48.365154 (50 iterations in 0.01 seconds)
    # Iteration 200: error is 50.364318 (50 iterations in 0.01 seconds)
    # Iteration 250: error is 50.683345 (50 iterations in 0.01 seconds)
    # Iteration 300: error is 1.136895 (50 iterations in 0.01 seconds)
    # Iteration 350: error is 0.733483 (50 iterations in 0.01 seconds)
    # Iteration 400: error is 0.625430 (50 iterations in 0.01 seconds)
    # Iteration 450: error is 0.625001 (50 iterations in 0.01 seconds)
    # Iteration 500: error is 0.624999 (50 iterations in 0.01 seconds)
    # Iteration 550: error is 0.624999 (50 iterations in 0.01 seconds)
    # Iteration 600: error is 0.624999 (50 iterations in 0.01 seconds)
    # Iteration 650: error is 0.624986 (50 iterations in 0.01 seconds)
    # Iteration 700: error is 0.624986 (50 iterations in 0.01 seconds)
    # Iteration 750: error is 0.625000 (50 iterations in 0.01 seconds)
    # Iteration 800: error is 0.624986 (50 iterations in 0.01 seconds)
    # Iteration 850: error is 0.624986 (50 iterations in 0.01 seconds)
    # Iteration 900: error is 0.625000 (50 iterations in 0.01 seconds)
    # Iteration 950: error is 0.625000 (50 iterations in 0.01 seconds)
    # Iteration 1000: error is 0.624999 (50 iterations in 0.01 seconds)
    # Fitting performed in 0.17 seconds.
    # Read the 115 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 37.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.01 seconds (sparsity = 0.989792)!
    # Learning embedding...
    # Iteration 50: error is 49.065562 (50 iterations in 0.01 seconds)
    # Iteration 100: error is 49.210910 (50 iterations in 0.01 seconds)
    # Iteration 150: error is 49.974686 (50 iterations in 0.01 seconds)
    # Iteration 200: error is 51.228638 (50 iterations in 0.01 seconds)
    # Iteration 250: error is 49.220168 (50 iterations in 0.01 seconds)
    # Iteration 300: error is 1.749707 (50 iterations in 0.01 seconds)
    # Iteration 350: error is 1.094096 (50 iterations in 0.01 seconds)
    # Iteration 400: error is 0.868036 (50 iterations in 0.01 seconds)
    # Iteration 450: error is 0.845715 (50 iterations in 0.01 seconds)
    # Iteration 500: error is 0.845736 (50 iterations in 0.01 seconds)
    # Iteration 550: error is 0.845737 (50 iterations in 0.01 seconds)
    # Iteration 600: error is 0.845737 (50 iterations in 0.01 seconds)
    # Iteration 650: error is 0.845737 (50 iterations in 0.01 seconds)
    # Iteration 700: error is 0.845737 (50 iterations in 0.01 seconds)
    # Iteration 750: error is 0.845736 (50 iterations in 0.01 seconds)
    # Iteration 800: error is 0.845736 (50 iterations in 0.01 seconds)
    # Iteration 850: error is 0.845737 (50 iterations in 0.01 seconds)
    # Iteration 900: error is 0.845737 (50 iterations in 0.01 seconds)
    # Iteration 950: error is 0.845737 (50 iterations in 0.01 seconds)
    # Iteration 1000: error is 0.845736 (50 iterations in 0.01 seconds)
    # Fitting performed in 0.16 seconds.
    # Read the 24 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 6.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.00 seconds (sparsity = 0.881944)!
    # Learning embedding...
    # Iteration 50: error is 70.976567 (50 iterations in 0.00 seconds)
    # Iteration 100: error is 72.024095 (50 iterations in 0.00 seconds)
    # Iteration 150: error is 72.360917 (50 iterations in 0.00 seconds)
    # Iteration 200: error is 71.631815 (50 iterations in 0.00 seconds)
    # Iteration 250: error is 73.191667 (50 iterations in 0.00 seconds)
    # Iteration 300: error is 6.587625 (50 iterations in 0.00 seconds)
    # Iteration 350: error is 3.246346 (50 iterations in 0.00 seconds)
    # Iteration 400: error is 2.966947 (50 iterations in 0.00 seconds)
    # Iteration 450: error is 2.960004 (50 iterations in 0.00 seconds)
    # Iteration 500: error is 2.959519 (50 iterations in 0.00 seconds)
    # Iteration 550: error is 2.959505 (50 iterations in 0.00 seconds)
    # Iteration 600: error is 2.959505 (50 iterations in 0.00 seconds)
    # Iteration 650: error is 2.959505 (50 iterations in 0.00 seconds)
    # Iteration 700: error is 2.959505 (50 iterations in 0.00 seconds)
    # Iteration 750: error is 2.959504 (50 iterations in 0.00 seconds)
    # Iteration 800: error is 2.959503 (50 iterations in 0.00 seconds)
    # Iteration 850: error is 2.959503 (50 iterations in 0.00 seconds)
    # Iteration 900: error is 2.959502 (50 iterations in 0.00 seconds)
    # Iteration 950: error is 2.959501 (50 iterations in 0.00 seconds)
    # Iteration 1000: error is 2.959500 (50 iterations in 0.00 seconds)
    # Fitting performed in 0.02 seconds.
    # Read the 10 x 6 data matrix successfully!
    # OpenMP is working. 1 threads.
    # Using no_dims = 1, perplexity = 2.000000, and theta = 0.100000
    # Computing input similarities...
    # Building tree...
    # Done in 0.00 seconds (sparsity = 0.760000)!
    # Learning embedding...
    # Iteration 50: error is 91.006955 (50 iterations in 0.00 seconds)
    # Iteration 100: error is 69.804871 (50 iterations in 0.00 seconds)
    # Iteration 150: error is 82.678065 (50 iterations in 0.00 seconds)
    # Iteration 200: error is 82.559507 (50 iterations in 0.00 seconds)
    # Iteration 250: error is 130.669903 (50 iterations in 0.00 seconds)
    # Iteration 300: error is 2.313263 (50 iterations in 0.00 seconds)
    # Iteration 350: error is 2.268197 (50 iterations in 0.00 seconds)
    # Iteration 400: error is 2.263634 (50 iterations in 0.00 seconds)
    # Iteration 450: error is 2.262503 (50 iterations in 0.00 seconds)
    # Iteration 500: error is 2.262447 (50 iterations in 0.00 seconds)
    # Iteration 550: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 600: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 650: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 700: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 750: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 800: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 850: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 900: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 950: error is 2.262446 (50 iterations in 0.00 seconds)
    # Iteration 1000: error is 2.262446 (50 iterations in 0.00 seconds)
    # Fitting performed in 0.00 seconds.

This is the version of R and the packages that were used to generate
these results.

``` r
sessionInfo()
# R version 3.6.2 (2019-12-12)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] fastTopics_0.4-23       SeuratWrappers_0.3.2    pbmc3k.SeuratData_3.1.4
# [4] SeuratData_0.2.1        Seurat_3.2.3           
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_1.4-1      deldir_0.1-29        
#   [4] ggridges_0.5.2        spatstat.data_1.4-3   farver_2.0.1         
#   [7] leiden_0.3.3          listenv_0.8.0         remotes_2.1.0        
#  [10] MatrixModels_0.4-1    ggrepel_0.9.0         fansi_0.4.0          
#  [13] codetools_0.2-16      splines_3.6.2         knitr_1.26           
#  [16] polyclip_1.10-0       zeallot_0.1.0         jsonlite_1.6         
#  [19] mcmc_0.9-6            ica_1.0-2             cluster_2.1.0        
#  [22] png_0.1-7             uwot_0.1.10           shiny_1.4.0          
#  [25] sctransform_0.3.2     BiocManager_1.30.10   compiler_3.6.2       
#  [28] httr_1.4.2            backports_1.1.5       assertthat_0.2.1     
#  [31] Matrix_1.2-18         fastmap_1.0.1         lazyeval_0.2.2       
#  [34] cli_2.0.0             later_1.0.0           prettyunits_1.1.1    
#  [37] htmltools_0.4.0       quantreg_5.54         tools_3.6.2          
#  [40] rsvd_1.0.2            igraph_1.2.5          coda_0.19-3          
#  [43] gtable_0.3.0          glue_1.3.1            RANN_2.6.1           
#  [46] reshape2_1.4.3        dplyr_0.8.3           rappdirs_0.3.1       
#  [49] Rcpp_1.0.5            spatstat_1.64-1       scattermore_0.7      
#  [52] vctrs_0.2.1           nlme_3.1-142          lmtest_0.9-38        
#  [55] xfun_0.11             stringr_1.4.0         globals_0.13.0       
#  [58] mime_0.8              miniUI_0.1.1.1        lifecycle_0.1.0      
#  [61] irlba_2.3.3           goftest_1.2-2         future_1.18.0        
#  [64] MASS_7.3-51.4         zoo_1.8-7             scales_1.1.0         
#  [67] hms_0.5.2             promises_1.1.0        spatstat.utils_1.17-0
#  [70] parallel_3.6.2        SparseM_1.78          RColorBrewer_1.1-2   
#  [73] yaml_2.2.0            reticulate_1.16       pbapply_1.4-3        
#  [76] gridExtra_2.3         ggplot2_3.3.0         rpart_4.1-15         
#  [79] stringi_1.4.3         rlang_0.4.5           pkgconfig_2.0.3      
#  [82] matrixStats_0.56.0    evaluate_0.14         lattice_0.20-38      
#  [85] ROCR_1.0-11           purrr_0.3.3           tensor_1.5           
#  [88] labeling_0.3          patchwork_1.0.1       htmlwidgets_1.5.1    
#  [91] cowplot_1.0.0         tidyselect_0.2.5      RcppAnnoy_0.0.18     
#  [94] plyr_1.8.5            magrittr_1.5          R6_2.4.1             
#  [97] pillar_1.4.3          mgcv_1.8-31           fitdistrplus_1.1-1   
# [100] survival_3.1-8        abind_1.4-5           tibble_2.1.3         
# [103] future.apply_1.6.0    crayon_1.3.4          KernSmooth_2.23-16   
# [106] plotly_4.9.2          rmarkdown_2.3         progress_1.2.2       
# [109] grid_3.6.2            data.table_1.12.8     digest_0.6.23        
# [112] xtable_1.8-4          tidyr_1.0.0           httpuv_1.5.2         
# [115] MCMCpack_1.4-5        RcppParallel_4.4.2    munsell_0.5.0        
# [118] viridisLite_0.3.0     quadprog_1.5-8
```
