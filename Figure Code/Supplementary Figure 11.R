library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(scran)
library(scRNAseq)
library(BBmisc)

tabulafacsabs <- readRDS("Tabula FACS/tabulafacsabs.rds")
tabulafacsabs@meta.data$logRNA <- log2(tabulafacsabs@meta.data$nCount_RNA + 1)

tabulafacsabs <- NormalizeData(tabulafacsabs, normalization.method = "RC")
tabulafacsabs@assays$RNA@data <- as.matrix(log2(tabulafacsabs@assays$RNA@counts + 1))
tabulafacsabs <- FindVariableFeatures(tabulafacsabs, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabulafacsabs)
tabulafacsabs <- ScaleData(tabulafacsabs, features = all.genes)
tabulafacsabs <- RunPCA(tabulafacsabs, features = VariableFeatures(object = tabulafacsabs))
tabulafacsabs <- FindNeighbors(tabulafacsabs, dims = 1:15)
tabulafacsabs <- FindClusters(tabulafacsabs, resolution = 0.4)
tabulafacsabs <- RunUMAP(tabulafacsabs, dims = 1:15, min.dist = 0.75)

source("https://raw.githubusercontent.com/jumphone/Vector/master/Vector.R")
VEC <- tabulafacsabs@reductions$umap@cell.embeddings
rownames(VEC) <- colnames(tabulafacsabs)
PCA <- tabulafacsabs@reductions$pca@cell.embeddings
PCA <- vector.rankPCA(PCA)

OUT <- vector.buildGrid(VEC, N = 150, SHOW = TRUE)
OUT <- vector.buildNet(OUT, CUT = 1, SHOW = TRUE)
OUT <- vector.getValue(OUT, PCA, SHOW = TRUE)
tabulafacsabs@meta.data$QP <- OUT$VALUE

library(BBmisc)
cols <- read.csv(file = "tabulafacscols.csv")

agg <- as.data.frame(tabulafacsabs@meta.data$logRNA)
colnames(agg) <- "logRNA"
agg[, "Celltype"] <- tabulafacsabs@meta.data$celltype
agg[, "features"] <- tabulafacsabs@meta.data$nFeature_RNA
agg[, "RNA"] <- tabulafacsabs@meta.data$nCount_RNA
agg[, "QP"] <- tabulafacsabs@meta.data$QP
agg[, "Organ"] <- tabulafacsabs@meta.data$tissue
agg[, "Organcell"] <- paste(agg$Celltype, "/", agg$Organ, sep = "")
agg$Celltype <- NULL
agg$Organ <- NULL

agg <- aggregate(. ~ Organcell, agg, mean)
rownames(agg) <- agg$Organcell
agg[, "Organcell"] <- NULL
agg <- agg[-c(1:3), ]

# Coefficient of Variation
co.var <- function(x) {
  (
    100 * sd(x) / mean(x)
  )
}

agg_CV <- as.data.frame(tabulafacsabs@meta.data$nCount_RNA)
colnames(agg_CV) <- "RNA"
agg_CV[, "Celltype"] <- tabulafacsabs@meta.data$celltype
agg_CV[, "Organ"] <- tabulafacsabs@meta.data$tissue
agg_CV[, "Organcell"] <- paste(agg_CV$Celltype, "/", agg_CV$Organ, sep = "")
agg_CV$Celltype <- NULL
agg_CV$Organ <- NULL

agg_CV <- as.data.frame(aggregate(. ~ Organcell, agg_CV, function(x) co.var(x)))
rownames(agg_CV) <- agg_CV$Organcell
agg_CV[, "Organcell"] <- NULL
agg_CV <- as.data.frame(agg_CV[-c(1:3), ])
colnames(agg_CV) <- "CV"

agg <- cbind(agg, agg_CV)
agg[, "organcell"] <- rownames(agg)
agg[, "Organ"] <- sapply(X = strsplit(agg$organcell, split = "/"), FUN = "[", 2)
agg[, "Celltype"] <- sapply(X = strsplit(agg$organcell, split = "/"), FUN = "[", 1)

agg[, "colors"] <- cols$colors[match(agg$Organ, cols$Organs)]
agg$organcell <- paste(agg$Celltype, " (", agg$Organ, ")", sep = "")
rownames(agg) <- agg$organcell

ggplot(agg, aes(x = CV, y = features)) +
  geom_point(aes(color = colors), size = 2) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(x = "CV (Transcripts)", y = "Feature Counts") +
  stat_cor(method = "spearman", label.x.npc = "centre", label.y.npc = "bottom") +
  NoLegend()

ggplot(agg, aes(x = CV, y = QP)) +
  geom_point(aes(color = colors), size = 2) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(x = "CV (Transcripts)", y = "QP Score") +
  stat_cor(method = "spearman", label.x.npc = "centre", label.y.npc = "bottom") +
  NoLegend()

ggplot(agg, aes(x = CV, y = logRNA)) +
  geom_point(aes(color = colors), size = 2) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(x = "CV (Transcripts)", y = "log2 Transcripts") +
  stat_cor(method = "spearman", label.x.npc = "centre", label.y.npc = "bottom") +
  NoLegend()
