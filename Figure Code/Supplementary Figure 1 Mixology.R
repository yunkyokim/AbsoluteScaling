library(Seurat)
library(ggplot2)
library(BBmisc)
library(data.table)
library(dplyr)
library(scran)
library(scater)
library(tidyverse)

# Panels A-D
# Data import
scattercol <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1")

# For Sort-Seq
umi_raw <- read.csv(file = paste0("Mixology Sort-seq/UMI Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
ercc_raw <- read.csv(file = paste0("Mixology Sort-seq/ERCC Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
barcodes <- read.csv(file = paste0("Mixology Sort-seq/barcodes.csv"), sep = ",", header = TRUE, row.names = 1)
# For CEL-Seq
umi_raw <- read.csv(file = paste0("Mixology CEL-seq/UMI Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
ercc_raw <- read.csv(file = paste0("Mixology CEL-seq/ERCC Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
barcodes <- read.csv(file = paste0("Mixology CEL-seq/barcodes.csv"), sep = ",", header = TRUE, row.names = 1)

# ERCC Object
ERCC_Abs2 <- CreateSeuratObject(counts = ercc_raw)
ERCC_Abs2 <- AddMetaData(ERCC_Abs2, metadata = barcodes$mRNA_amount, col.name = "mRNAConcentration")
ERCC_Abs2 <- AddMetaData(ERCC_Abs2, PercentageFeatureSet(ERCC_Abs2, pattern = "^ERCC-"), col.name = "ERCC")
ERCC_Abs2 <- subset(ERCC_Abs2, subset = nCount_RNA > 200 & ERCC > 0.1)

filteredcounts <- GetAssayData(object = ERCC_Abs2, assay = "RNA", slot = "data")
scran.data <- SingleCellExperiment(assays = list(counts = filteredcounts))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <- splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))

scran.data <- computeSpikeFactors(scran.data, "ERCC")
scran.data <- logNormCounts(scran.data, log = FALSE)

ERCC_Abs <- CreateSeuratObject(counts = as.sparse(x = assay(scran.data, "normcounts")))
ERCC_Abs@meta.data$mRNAConcentration <- ERCC_Abs2@meta.data$mRNAConcentration
ERCC_Abs <- AddMetaData(ERCC_Abs, PercentageFeatureSet(ERCC_Abs, pattern = "^ERCC-"), col.name = "ERCC")
rm(ERCC_Abs2)

ERCC_Abs@assays$RNA@data <- as.matrix(log2(ERCC_Abs@assays$RNA@counts + 1))
Idents(object = ERCC_Abs) <- ERCC_Abs@meta.data$mRNAConcentration

# UMI Object
UMI_Abs <- CreateSeuratObject(counts = umi_raw)
UMI_Abs <- AddMetaData(UMI_Abs, metadata = barcodes$mRNA_amount, col.name = "mRNAConcentration")
UMI_Abs <- subset(UMI_Abs, subset = nCount_RNA > 200)

UMI_Abs@assays$RNA@data <- as.matrix(log2(UMI_Abs@assays$RNA@counts + 1))
Idents(object = UMI_Abs) <- UMI_Abs@meta.data$mRNAConcentration

# Comparison Scatter
UMI_Abs_transcripts <- data.frame("UMI Transcripts" = UMI_Abs@meta.data$nCount_RNA, "Names" = colnames(UMI_Abs@assays$RNA@data), "Value" = UMI_Abs@meta.data$mRNAConcentration)
ERCC_Abs_transcripts <- data.frame("ERCC Transcripts" = ERCC_Abs@meta.data$nCount_RNA, "Names" = colnames(ERCC_Abs@assays$RNA@data), "Value" = ERCC_Abs@meta.data$mRNAConcentration)
comp_Abs <- merge(UMI_Abs_transcripts, ERCC_Abs_transcripts, by = "Names")
comp_Abs$UMI.Transcripts <- normalize(comp_Abs$UMI.Transcripts, method = "range", range = c(0, 10))
comp_Abs$ERCC.Transcripts <- normalize(comp_Abs$ERCC.Transcripts, method = "range", range = c(0, 10))

library(ggpubr)
ggplot(comp_Abs, aes(x = ERCC.Transcripts, y = UMI.Transcripts)) +
  geom_point(size = 1.5, aes(color = as.factor(Value.x))) +
  geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "UMI-Normalization", x = "ERCC-Normalization") +
  scale_color_manual(values = c("seagreen4", "seagreen3", "seagreen2", "seagreen1")) +
  xlim(0, 10) +
  ylim(0, 10) +
  stat_regline_equation(size = 4.5) +
  NoLegend()


# Subsampling Analysis
# ERCC Subsampling
ercc1 <- CreateSeuratObject(counts = ercc_raw)
ercc1 <- AddMetaData(ercc1, metadata = barcodes$mRNA_amount, col.name = "mRNAConcentration")
ercc1 <- AddMetaData(ercc1, PercentageFeatureSet(ercc1, pattern = "^ERCC-"), col.name = "ERCC")
ercc1 <- subset(ercc1, subset = nCount_RNA > 200 & ERCC > 0.1)

ercc2 <- subset(ercc1, downsample = 143)
ercc3 <- subset(ercc1, downsample = 95)
ercc4 <- subset(ercc1, downsample = 47)
ercc5 <- subset(ercc1, downsample = 19)

auto_ercc <- function(ercc_dataset) {
  filteredcounts <- GetAssayData(object = ercc_dataset, assay = "RNA", slot = "data")
  scran.data <- SingleCellExperiment(assays = list(counts = filteredcounts))
  is.spike <- grepl("^ERCC-", rownames(scran.data))
  scran.data <- splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))

  scran.data <- computeSpikeFactors(scran.data, "ERCC")
  scran.data <- logNormCounts(scran.data, log = FALSE)

  ERCC_Abs <- CreateSeuratObject(counts = as.sparse(x = assay(scran.data, "normcounts")))
  ERCC_Abs@meta.data$mRNAConcentration <- ercc_dataset@meta.data$mRNAConcentration
  ERCC_Abs <- AddMetaData(ERCC_Abs, PercentageFeatureSet(ERCC_Abs, pattern = "^ERCC-"), col.name = "ERCC")

  ERCC_Abs@assays$RNA@data <- as.matrix(log2(ERCC_Abs@assays$RNA@counts + 1))
  Idents(object = ERCC_Abs) <- ERCC_Abs@meta.data$mRNAConcentration
  return(ERCC_Abs)
}

ercc1 <- auto_ercc(ercc1)
scattercol <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1")
FeatureScatter(ercc1, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

ercc2 <- auto_ercc(ercc2)
scattercol <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1")
FeatureScatter(ercc2, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

ercc3 <- auto_ercc(ercc3)
scattercol <- c("seagreen4", "seagreen2", "seagreen3", "seagreen1")
FeatureScatter(ercc3, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

ercc4 <- auto_ercc(ercc4)
scattercol <- c("seagreen2", "seagreen3", "seagreen1", "seagreen4")
FeatureScatter(ercc4, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

ercc5 <- auto_ercc(ercc5)
scattercol <- c("seagreen2", "seagreen3", "seagreen1", "seagreen4")
FeatureScatter(ercc5, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()


# UMI Subsampling
umi1 <- CreateSeuratObject(counts = umi_raw)
umi1 <- AddMetaData(umi1, metadata = barcodes$mRNA_amount, col.name = "mRNAConcentration")
umi1 <- subset(umi1, subset = nCount_RNA > 200)

umi2 <- subset(umi1, downsample = 287)
umi3 <- subset(umi1, downsample = 191)
umi4 <- subset(umi1, downsample = 95)
umi5 <- subset(umi1, downsample = 38)

umi1@assays$RNA@data <- as.matrix(log2(umi1@assays$RNA@counts + 1))
Idents(object = umi1) <- umi1@meta.data$mRNAConcentration
scattercol <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1")
FeatureScatter(umi1, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

umi2@assays$RNA@data <- as.matrix(log2(umi2@assays$RNA@counts + 1))
Idents(object = umi2) <- umi2@meta.data$mRNAConcentration
scattercol <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1")
FeatureScatter(umi2, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

umi3@assays$RNA@data <- as.matrix(log2(umi3@assays$RNA@counts + 1))
Idents(object = umi3) <- umi3@meta.data$mRNAConcentration
scattercol <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1")
FeatureScatter(umi3, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

umi4@assays$RNA@data <- as.matrix(log2(umi4@assays$RNA@counts + 1))
Idents(object = umi4) <- umi4@meta.data$mRNAConcentration
scattercol <- c("seagreen2", "seagreen3", "seagreen1", "seagreen4")
FeatureScatter(umi4, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()

umi5@assays$RNA@data <- as.matrix(log2(umi5@assays$RNA@counts + 1))
Idents(object = umi5) <- umi5@meta.data$mRNAConcentration
scattercol <- c("seagreen2", "seagreen3", "seagreen1", "seagreen4")
FeatureScatter(umi5, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  NoLegend()
