library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(scran)
library(scRNAseq)
library(scater)
library(BBmisc)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(viridis)
library(DoubletFinder)

# Dataset Prep
gonadE11 <- Read10X(data.dir = "Niu2020/E11.5/", gene.column = 2, unique.features = TRUE)
gonad_E11 <- CreateSeuratObject(gonadE11)
gonadE12 <- Read10X(data.dir = "Niu2020/E12.5/", gene.column = 2, unique.features = TRUE)
gonad_E12 <- CreateSeuratObject(gonadE12)
gonadE14 <- Read10X(data.dir = "Niu2020/E14.5/", gene.column = 2, unique.features = TRUE)
gonad_E14 <- CreateSeuratObject(gonadE14)
gonadE16 <- Read10X(data.dir = "Niu2020/E16.5/", gene.column = 2, unique.features = TRUE)
gonad_E16 <- CreateSeuratObject(gonadE16)
gonadE18 <- Read10X(data.dir = "Niu2020//E18.5/", gene.column = 2, unique.features = TRUE)
gonad_E18 <- CreateSeuratObject(gonadE18)
gonadP1 <- Read10X(data.dir = "Niu2020/P1/", gene.column = 2, unique.features = TRUE)
gonad_P1 <- CreateSeuratObject(gonadP1)
gonadP5 <- Read10X(data.dir = "Niu2020/P5/", gene.column = 2, unique.features = TRUE)
gonad_P5 <- CreateSeuratObject(gonadP5)
rm(gonadE11, gonadE12, gonadE14, gonadE16, gonadE18, gonadP1, gonadP5)

# Bulk Doublet Detection
timepoints <- as.list(ls())

# Perform Doublet Detection
auto_doublet_detection <- function(seurat_dataset, est_doublet) {
  # Log - Normalized Doublet Detection
  seurat_dataset <- subset(seurat_dataset, subset = nFeature_RNA > 1000 & nCount_RNA > 500 & nCount_RNA < 60000)
  seurat_dataset <- NormalizeData(seurat_dataset, normalization.method = "LogNormalize")
  seurat_dataset < -ScaleData(seurat_dataset)
  seurat_dataset < -FindVariableFeatures(seurat_dataset, selection.method = "vst", nfeatures = 2000)
  seurat_dataset < -RunPCA(seurat_dataset)
  seurat_dataset < -FindNeighbors(seurat_dataset, dims = 1:10)
  seurat_dataset < -FindClusters(seurat_dataset, resolution = 0.4)
  seurat_dataset < -RunUMAP(seurat_dataset, dims = 1:10)

  # Doublet Detection
  sweep.res < -paramSweep_v3(seurat_dataset, PCs = 1:10, sct = FALSE)
  sweep.stats < -summarizeSweep(sweep.res, GT = FALSE)
  log.bcmvn < -(as.data.frame(find.pK(sweep.stats)))
  log.bcmvn$pK < -as.numeric(as.character(log.bcmvn$pK))
  max_pk <- log.bcmvn$pK[[which(grepl(max(data.matrix(log.bcmvn$BCmetric)), log.bcmvn$BCmetric))]]

  annotations <- seurat_dataset @meta.data$seurat_clusters
  homotypic.prop < -modelHomotypic(annotations)
  est_doublet <- 0.04
  nExp_poi < -round(est_doublet * length(seurat_dataset$orig.ident))
  nExp_poi.adj < -round(nExp_poi * (1 - homotypic.prop))
  seurat_dataset < -doubletFinder_v3(seurat_dataset, PCs = 1:10, pN = 0.25, pK = max_pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  meta_doublets <- paste("DF.classifications_0.25_", as.character(max_pk), "_", as.character(nExp_poi.adj), sep = "")

  # Subset
  seurat_dataset <- eval(parse(text = paste("subset(seurat_dataset, subset = ", meta_doublets, "== 'Singlet')", sep = "")))
  rm(log.bcmvn, sweep.res, sweep.stats, annotations, homotypic.prop, nExp_poi, nExp_poi.adj, max_pk, est_doublet, meta_doublets)
  return(seurat_dataset)
}
for (organ_dataset in timepoints) {
  eval(parse(text = paste(organ_dataset, "= auto_doublet_detection(", organ_dataset, ", 0.04)", sep = "")))
}

# Merge
gonads <- merge(gonad_E11, y = c(gonad_E12, gonad_E14, gonad_E16, gonad_E18, gonad_P1, gonad_P5), add.cell.ids = c("E11.5", "E12.5", "E14.5", "E16.5", "E18.5", "P1", "P5"), project = "gonads")
rm(list = setdiff(ls(), "gonads"))

timepoint <- sapply(X = strsplit(colnames(gonads), split = "_"), FUN = "[", 1)
gonads <- AddMetaData(object = gonads, metadata = timepoint, col.name = "timepoint")
saveRDS(gonads, file = "Niu2020/gonadstotal.rds")

# Analysis
gonads <- readRDS("Niu2020/gonadstotal.rds")
gonads @meta.data$noclusters <- 0
gonads <- subset(gonads, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & nCount_RNA < 50000) # basic filtering
gonads <- subset(gonads, subset = Car2 == 0 & Cldn5 == 0) # remove immune and other cells

gonads @assays$RNA @data <- as.matrix(log2(gonads @assays$RNA @counts + 1))
gonads < -FindVariableFeatures(gonads)
all.genes < -rownames(gonads)
gonads < -ScaleData(gonads, features = all.genes)
gonads < -RunPCA(gonads, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
gonads < -FindNeighbors(gonads, dims = 1:15)
gonads < -FindClusters(gonads, resolution = 2)
gonads < -RunUMAP(gonads, dims = 1:15, min.dist = 0.75)

Oct4 <- as.data.frame(gonads @assays$RNA @counts["Pou5f1", ])
colnames(Oct4) <- "Oct4"
for (row in 1:nrow(Oct4)) {
  expr < -Oct4[row, "Oct4"]
  if (expr > 0.5) {
    Oct4[row, "status"] <- "Oct4+"
  } else {
    Oct4[row, "status"] <- "Oct4-"
  }
}
# assign PGCs by Oct4 expression
gonads @meta.data$germ <- Oct4$status

# Figure Panels
dimcols <- c("grey80", "violetred2")
dimcols2 <- c("grey80", "violetred2", "grey80", "violetred2", "grey80", "violetred2", "grey80", "violetred2", "grey80", "violetred2", "grey80", "violetred2", "grey80", "violetred2")
VlnPlot(gonads, features = c("nCount_RNA"), split.plot = TRUE, group.by = "timepoint", cols = dimcols, split.by = "germ", pt.size = 0) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  stat_summary(fun.y = median, geom = "point", size = 10, colour = dimcols2, shape = 95) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

gonads_log <- NormalizeData(gonads, normalization.method = "RC")
gonads_log <- CreateSeuratObject(counts = gonads_log @assays$RNA @data)
gonads_log @assays$RNA @data <- as.matrix(log2(gonads_log @assays$RNA @counts + 1))
gonads_log @meta.data$nCount_logRNA <- colSums(gonads_log @assays$RNA @data)
gonads_log @meta.data$germ <- gonads @meta.data$germ
gonads_log @meta.data$timepoint <- gonads @meta.data$timepoint

VlnPlot(gonads_log, features = c("nCount_logRNA"), split.plot = TRUE, group.by = "timepoint", cols = dimcols, split.by = "germ", pt.size = 0) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  stat_summary(fun.y = median, geom = "point", size = 10, colour = dimcols2, shape = 95) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

saveRDS(gonads_log, "Niu2020/gonads_glo.rds")
saveRDS(gonads, "Niu2020/gonads_abs.rds")

# Cell Cycle Scoring
library(stringr)
s.genes < -str_to_title(cc.genes$s.genes)
g2m.genes < -str_to_title(cc.genes$g2m.genes)
gonads_log < -CellCycleScoring(gonads_log, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

VlnPlot(gonads_log, features = c("S.Score"), split.plot = TRUE, group.by = "timepoint", cols = dimcols, split.by = "germ", pt.size = 0) +
  stat_summary(fun.y = median, geom = "point", size = 10, colour = dimcols2, shape = 95) +
  ylab("S-Phase Score") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

VlnPlot(gonads_log, features = c("G2M.Score"), split.plot = TRUE, group.by = "timepoint", cols = dimcols, split.by = "germ", pt.size = 0) +
  stat_summary(fun.y = median, geom = "point", size = 10, colour = dimcols2, shape = 95) +
  ylab("G2/M-Phase Score") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# E12 .5 Panels
gonad_E12 <- subset(gonads, subset = timepoint == "E12.5")
gonad_E12 < -FindVariableFeatures(gonad_E12)
all.genes < -rownames(gonad_E12)
gonad_E12 < -ScaleData(gonad_E12, features = all.genes)
gonad_E12 < -RunPCA(gonad_E12, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
gonad_E12 < -FindNeighbors(gonad_E12, dims = 1:15)
gonad_E12 < -FindClusters(gonad_E12, resolution = 0.5)
gonad_E12 < -RunUMAP(gonad_E12, dims = 1:15, min.dist = 0.75)

Oct4 <- as.data.frame(gonad_E12 @assays$RNA @counts["Pou5f1", ])
colnames(Oct4) <- "Oct4"
for (row in 1:nrow(Oct4)) {
  expr < -Oct4[row, "Oct4"]

  if (expr > 0.3) {
    Oct4[row, "status"] <- "Oct4+"
  } else {
    Oct4[row, "status"] <- "Oct4-"
  }
}

gonad_E12 @meta.data$germ <- Oct4$status
gonad_E12_log <- NormalizeData(gonad_E12, normalization.method = "RC")
gonad_E12_log <- CreateSeuratObject(counts = gonad_E12_log @assays$RNA @data)
gonad_E12_log @assays$RNA @data <- as.matrix(log2(gonad_E12_log @assays$RNA @counts + 1))
gonad_E12_log @meta.data$germ <- gonad_E12 @meta.data$germ
gonad_E12_log < -FindVariableFeatures(gonad_E12_log)
all.genes < -rownames(gonad_E12_log)
gonad_E12_log < -ScaleData(gonad_E12_log, features = all.genes)
gonad_E12_log < -RunPCA(gonad_E12_log, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

gonad_E12_log < -FindNeighbors(gonad_E12_log, dims = 1:15)
gonad_E12_log < -FindClusters(gonad_E12_log, resolution = 0.5)
gonad_E12_log < -RunUMAP(gonad_E12_log, dims = 1:15, min.dist = 0.75)
gonad_E12_log @meta.data$nCount_logRNA <- colSums(gonad_E12_log @assays$RNA @data)

# Figure Panels
DimPlot(gonad_E12_log, reduction = "umap", group.by = "germ", cols = dimcols, pt.size = 0.1) +
  theme(aspect.ratio = 1, text = element_text(size = 14), plot.title = element_blank(), axis.text = element_text(size = 12), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + NoLegend()

DimPlot(gonad_E12, reduction = "umap", group.by = "germ", cols = dimcols, pt.size = 0.1) +
  theme(aspect.ratio = 1, text = element_text(size = 14), plot.title = element_blank(), axis.text = element_text(size = 12), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + NoLegend()

FeaturePlot(gonad_E12_log, features = c("nCount_logRNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

FeaturePlot(gonad_E12, features = c("nCount_RNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )
