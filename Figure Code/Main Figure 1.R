library(Seurat)
library(ggplot2)
library(BBmisc)
library(data.table)
library(dplyr)
library(scran)
library(scater)
library(tidyverse)

# Panel C
# 10 Mixture

outersect <-
  function(x, y) {
    sort(c(
      setdiff(x, y),
      setdiff(y, x)
    ))
  }

# Dataset Processing
hg.mix <- Read10X(data.dir = "10X Mixture/hg19v2/")
mm.mix <- Read10X(data.dir = "10X Mixture/mm10v2/")
joint.bcs <- intersect(colnames(hg.mix), colnames(mm.mix))
hg.bcs <- intersect(outersect(colnames(hg.mix), colnames(mm.mix)), colnames(hg.mix))
mm.bcs <- intersect(outersect(colnames(hg.mix), colnames(mm.mix)), colnames(mm.mix))

hg.joint <- hg.mix[, joint.bcs]
mm.joint <- mm.mix[, joint.bcs]
hgmm.joint <- rbind(hg.joint, mm.joint)

hg.unique <- hg.mix[, hg.bcs]
mm.unique <- mm.mix[, mm.bcs]

hg.seurat <- CreateSeuratObject(counts = hg.unique, project = "Human")
mm.seurat <- CreateSeuratObject(counts = mm.unique, project = "Mouse")
hgmm.seurat <- CreateSeuratObject(counts = hgmm.joint, project = "Mixed")

mixture <- merge(hgmm.seurat, y = c(hg.seurat, mm.seurat), add.cell.ids = c("Mixed", "Human", "Mouse"), project = "Mixed Doublets")
mixture @meta.data$noclusters <- 0
mixture <- subset(mixture, subset = nFeature_RNA > 500 & nCount_RNA > 1000)
saveRDS(mixture, "Mixture/mixture.subset")

# Plot
vlncol <- c("hotpink4", "hotpink3", "hotpink2")
VlnPlot(mixture, features = c("nCount_RNA"), group.by = "orig.ident", cols = vlncol, ncol = 1, pt.size = 0, sort = TRUE) +
  stat_summary(fun.y = median, geom = "point", size = 15, colour = "grey27", shape = 95) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  ylab("Transcripts") +
  xlab("") +
  NoLegend() +
  ggtitle("")

# Cell Hashing
# Dataset Processing
hashing.umi <- readRDS(file = "Cell Hashing/pbmc_umi_mtx.rds")
hashing.hto <- readRDS(file = "Cell Hashing/pbmc_hto_mtx.rds")
joint.bcs <- intersect(colnames(hashing.umi), colnames(hashing.hto))
hashing.umi <- hashing.umi[, joint.bcs]
hashing.hto <- as.matrix(hashing.hto[, joint.bcs])

hashing <- CreateSeuratObject(counts = hashing.umi)
hashing[["HTO"]] <- CreateAssayObject(counts = hashing.hto)
hashing <- HTODemux(hashing, assay = "HTO", positive.quantile = 0.99)
table(hashing$HTO_classification.global)

Idents(hashing) <- "HTO_classification.global"
hashing.subset <- subset(hashing, idents = "Negative", invert = TRUE)
hashing.subset @meta.data$orig.ident <- hashing.subset @meta.data$HTO_classification.global
hashing.subset @meta.data$noclusters <- 0
hashing.subset <- subset(hashing.subset, subset = nFeature_RNA > 150 & nCount_RNA > 100 & nCount_RNA < 3000)
saveRDS(hashing.subset, "Cell Hashing/hashing.subset")

# Plot
vlncol <- c("hotpink4", "hotpink3")
VlnPlot(hashing.subset, features = c("nCount_RNA"), group.by = "orig.ident", cols = vlncol, ncol = 1, pt.size = 0) +
  stat_summary(fun.y = median, geom = "point", size = 15, colour = "grey27", shape = 95) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  ylab("Transcripts") +
  xlab("") +
  ggtitle("") +
  NoLegend()

# Demuxlet
# Dataset Processing
raw.data <- Read10X(data.dir = "Demuxlet/")
demuxlet <- CreateSeuratObject(counts = raw.data, project = "demuxlet")
cellinfo <- as.data.frame(read.table(file = "Demuxlet/demuxlet_calls.tsv", header = TRUE, row.names = 1))
demuxlet <- AddMetaData(demuxlet, cellinfo$Call, col.name = "orig.ident")
demuxlet @meta.data$noclusters <- 0

demuxlet <- subset(demuxlet, subset = nFeature_RNA > 200 & nCount_RNA > 100)
demuxlet <- subset(demuxlet, subset = orig.ident == "AMB", invert = TRUE)
cellstats <- demuxlet @meta.data

for (cellid in 1:nrow(cellstats)) {
  if (cellstats[cellid, "orig.ident"] == "SNG") {
    cellstats[cellid, "status"] <- "Singlet"
  }

  if (cellstats[cellid, "orig.ident"] == "DBL") {
    cellstats[cellid, "status"] <- "Doublet"
  }
}

demuxlet <- AddMetaData(demuxlet, cellstats$status, col.name = "status")
saveRDS(demuxlet, "Demuxlet/demuxlet.subset")

vlncol <- c("hotpink4", "hotpink3")
VlnPlot(demuxlet, features = c("nCount_RNA"), group.by = "status", cols = vlncol, ncol = 1, pt.size = 0) + ylab("Transcripts") + xlab("") + NoLegend() + ggtitle("") +
  stat_summary(fun.y = median, geom = "point", size = 15, colour = "grey27", shape = 95) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  NoLegend()

# Panels D - G
# Data import
scattercol <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1")

# For Sort - Seq
umi_raw <- read.csv(file = paste0("Mixology Sort-seq/UMI Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
ercc_raw <- read.csv(file = paste0("Mixology Sort-seq/ERCC Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
barcodes <- read.csv(file = paste0("Mixology Sort-seq/barcodes.csv"), sep = ",", header = TRUE, row.names = 1)
# For CEL - Seq
umi_raw <- read.csv(file = paste0("Mixology CEL-seq/UMI Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
ercc_raw <- read.csv(file = paste0("Mixology CEL-seq/ERCC Gene Count.csv"), sep = ",", header = TRUE, row.names = 1)
barcodes <- read.csv(file = paste0("Mixology CEL-seq/barcodes.csv"), sep = ",", header = TRUE, row.names = 1)

# UMI Global - scaling by relative counts
UMI_Glo2 <- CreateSeuratObject(counts = umi_raw)
UMI_Glo2 <- AddMetaData(UMI_Glo2, metadata = barcodes$mRNA_amount, col.name = "mRNAConcentration")
UMI_Glo2 <- subset(UMI_Glo2, subset = nCount_RNA > 200)

UMI_Glo2 <- NormalizeData(UMI_Glo2, normalization.method = "RC")
UMI_Glo2 @assays$RNA @data <- as.matrix(log2(UMI_Glo2 @assays$RNA @data + 1))
UMI_Glo <- CreateSeuratObject(counts = UMI_Glo2 @assays$RNA @data)
UMI_Glo @assays$RNA @counts <- UMI_Glo2 @assays$RNA @counts
UMI_Glo @assays$RNA @data <- UMI_Glo2 @assays$RNA @data
UMI_Glo @meta.data$mRNAConcentration <- UMI_Glo2 @meta.data$mRNAConcentration
rm(UMI_Glo2)

Idents(object = UMI_Glo) <- UMI_Glo @meta.data$mRNAConcentration

# Panel
FeatureScatter(UMI_Glo, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) # 600 x 350

# Transcriptome Curves
UMI_Glo_curve <- as.data.frame(t(as.data.frame(UMI_Glo @assays$RNA @data)))
UMI_Glo_curve[, "celltype"] <- UMI_Glo @active.ident
UMI_Glo_curve <- aggregate(. ~ celltype, UMI_Glo_curve, mean)
rownames(UMI_Glo_curve) <- UMI_Glo_curve$celltype
UMI_Glo_curve$celltype <- NULL
UMI_Glo_curve <- as.data.frame(t(UMI_Glo_curve))
UMI_Glo_curve[, "total"] <- rowSums(UMI_Glo_curve)
UMI_Glo_curve[, "genes"] <- rownames(UMI_Glo_curve)
UMI_Glo_curve <- UMI_Glo_curve[order(UMI_Glo_curve$total), ]

order.total <- order(UMI_Glo_curve$total, decreasing = TRUE)
UMI_Glo_curve$rank <- NA
UMI_Glo_curve$rank[order.total] <- 1:nrow(UMI_Glo_curve)

UMI_Glo_curve2 <- UMI_Glo_curve[UMI_Glo_curve$rank <= 20000, ]
UMI_Glo_tidy <- UMI_Glo_curve2 % > %
  select(`3.75`, `7.5`, `15`, `30`, rank) % > %
  gather(key = "variable", value = "value", -rank)

# Plot
ggplot(UMI_Glo_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), method = "gam", size = 1.5, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = c("seagreen2", "seagreen4", "seagreen1", "seagreen3")) +
  NoLegend()

# UMI Absolute Scaling by Raw UMIs
UMI_Abs <- CreateSeuratObject(counts = umi_raw)
UMI_Abs <- AddMetaData(UMI_Abs, metadata = barcodes$mRNA_amount, col.name = "mRNAConcentration")
UMI_Abs <- subset(UMI_Abs, subset = nCount_RNA > 200)

UMI_Abs @assays$RNA @data <- as.matrix(log2(UMI_Abs @assays$RNA @counts + 1))
Idents(object = UMI_Abs) <- UMI_Abs @meta.data$mRNAConcentration

# Plot
FeatureScatter(UMI_Abs, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )
xlab("Transcripts") +
  ylab("mRNA (pg)")

# Transcriptome Curves
UMI_Abs_curve <- as.data.frame(t(as.data.frame(UMI_Abs @assays$RNA @data)))
UMI_Abs_curve[, "celltype"] <- UMI_Abs @active.ident
UMI_Abs_curve <- aggregate(. ~ celltype, UMI_Abs_curve, mean)
rownames(UMI_Abs_curve) <- UMI_Abs_curve$celltype
UMI_Abs_curve$celltype <- NULL
UMI_Abs_curve <- as.data.frame(t(UMI_Abs_curve))
UMI_Abs_curve[, "total"] <- rowSums(UMI_Abs_curve)
UMI_Abs_curve[, "genes"] <- rownames(UMI_Abs_curve)
UMI_Abs_curve <- UMI_Abs_curve[order(UMI_Abs_curve$total), ]

order.total <- order(UMI_Abs_curve$total, decreasing = TRUE)
UMI_Abs_curve$rank <- NA
UMI_Abs_curve$rank[order.total] <- 1:nrow(UMI_Abs_curve)

UMI_Abs_curve2 <- UMI_Abs_curve
UMI_Glo_curve2 <- UMI_Glo_curve2[UMI_Glo_curve2$rank <= 20000, ]
UMI_Abs_curve2$`3.75` <- log2(UMI_Abs_curve2$`3.75` + 1)
UMI_Abs_curve2$`7.5` <- log2(UMI_Abs_curve2$`7.5` + 1)
UMI_Abs_curve2$`15` <- log2(UMI_Abs_curve2$`15` + 1)
UMI_Abs_curve2$`30` <- log2(UMI_Abs_curve2$`30` + 1)

UMI_Abs_tidy <- UMI_Abs_curve2 % > %
  select(`3.75`, `7.5`, `15`, `30`, rank) % > %
  gather(key = "variable", value = "value", -rank)

# Plot
ggplot(UMI_Abs_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.5, alpha = 1, se = FALSE) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = c("seagreen2", "seagreen4", "seagreen1", "seagreen3")) +
  NoLegend()

# ERCC Absolute Scaling by scran
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
ERCC_Abs @meta.data$mRNAConcentration <- ERCC_Abs2 @meta.data$mRNAConcentration
ERCC_Abs <- AddMetaData(ERCC_Abs, PercentageFeatureSet(ERCC_Abs, pattern = "^ERCC-"), col.name = "ERCC")
rm(ERCC_Abs2)

ERCC_Abs @assays$RNA @data <- as.matrix(log2(ERCC_Abs @assays$RNA @counts + 1))
Idents(object = ERCC_Abs) <- ERCC_Abs @meta.data$mRNAConcentration

FeatureScatter(ERCC_Abs, "nCount_RNA", "mRNAConcentration", slot = "data", pt.size = 1.5, cols = scattercol) + xlab("Transcripts") + ylab("mRNA (pg)") +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) # 600 x 350

# Transcriptome Curve
ERCC_Abs_curve <- as.data.frame(t(as.data.frame(ERCC_Abs @assays$RNA @data)))
ERCC_Abs_curve[, "celltype"] <- ERCC_Abs @active.ident
ERCC_Abs_curve <- aggregate(. ~ celltype, ERCC_Abs_curve, mean)
rownames(ERCC_Abs_curve) <- ERCC_Abs_curve$celltype
ERCC_Abs_curve$celltype <- NULL
ERCC_Abs_curve <- as.data.frame(t(ERCC_Abs_curve))
ERCC_Abs_curve[, "total"] <- rowSums(ERCC_Abs_curve)
ERCC_Abs_curve[, "genes"] <- rownames(ERCC_Abs_curve)
ERCC_Abs_curve <- ERCC_Abs_curve[order(ERCC_Abs_curve$total), ]

order.total <- order(ERCC_Abs_curve$total, decreasing = TRUE)
ERCC_Abs_curve$rank <- NA
ERCC_Abs_curve$rank[order.total] <- 1:nrow(ERCC_Abs_curve)

ERCC_Abs_curve2 <- ERCC_Abs_curve[ERCC_Abs_curve$rank <= 20000, ]
ERCC_Abs_tidy <- ERCC_Abs_curve2 % > %
  select(`3.75`, `7.5`, `15`, `30`, rank) % > %
  gather(key = "variable", value = "value", -rank)

# Plot
ggplot(ERCC_Abs_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.5, alpha = 1, se = FALSE) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = c("seagreen2", "seagreen4", "seagreen1", "seagreen3")) +
  NoLegend()
