library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)
library(viridis)
library(stats)
library(ggridges)

# Global UMAPs
# Tabula FACS, global scaling + absolute transcripts projected
tabulafacs <- readRDS("Tabula FACS/tabulafacsabs.rds")
tabula_unnorm <- readRDS("Tabula FACS/tabulaunnorm.rds")
tabulafacs@meta.data$logRNA = log2(tabulafacs@meta.data$nCount_RNA+1)
tabula_unnorm@meta.data$logRNA <- tabulafacs@meta.data$logRNA
tabula_unnorm@meta.data$facsRNA <- tabulafacs@meta.data$nCount_RNA

tabula_unnorm <- NormalizeData(tabula_unnorm, normalization.method = "RC")
tabula_unnorm@assays$RNA@data <- as.matrix(log2(tabula_unnorm@assays$RNA@data + 1))
tabula_unnorm <- FindVariableFeatures(tabula_unnorm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabula_unnorm)
tabula_unnorm <- ScaleData(tabula_unnorm, features = all.genes)
tabula_unnorm <- RunPCA(tabula_unnorm, features = VariableFeatures(object = tabula_unnorm))
tabula_unnorm <- FindNeighbors(tabula_unnorm, dims = 1:15)
tabula_unnorm <- FindClusters(tabula_unnorm, resolution = 0.4)
tabula_unnorm <- RunUMAP(tabula_unnorm, dims = 1:15, min.dist = 0.75)

FeaturePlot(tabula_unnorm, features = c("logRNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (Smart-seq2)") # 450 x 350

DimPlot(tabula_unnorm, reduction = "umap", group.by = "tissue", pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (Smart-seq2)") # 450 x 350

# Tabula FACS, absolute scaling + absolute transcripts projected
tabulafacs <- readRDS("Tabula FACS/tabulafacsabs.rds")
tabulafacs@assays$RNA@data <- as.matrix(log2(tabulafacs@assays$RNA@counts + 1))
tabulafacs <- FindVariableFeatures(tabulafacs, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabulafacs)
tabulafacs <- ScaleData(tabulafacs, features = all.genes)
tabulafacs <- RunPCA(tabulafacs, features = VariableFeatures(object = tabulafacs))
tabulafacs <- FindNeighbors(tabulafacs, dims = 1:15)
tabulafacs <- FindClusters(tabulafacs, resolution = 0.4)
tabulafacs <- RunUMAP(tabulafacs, dims = 1:15, min.dist = 0.75)

FeaturePlot(tabulafacs, features = c("logRNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (Smart-seq2)") # 450 x 350

DimPlot(tabulafacs, reduction = "umap", group.by = "tissue", pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (Smart-seq2)") # 450 x 350

# Tabula 10X, global scaling + absolute transcripts projected
tabula10X <- readRDS("Tabula 10X/tabula10X.rds")
tabula10X@meta.data$logRNA <- log2(tabula10X@meta.data$nCount_RNA + 1)

tabula10X <- NormalizeData(tabula10X, normalization.method = "RC")
tabula10X@assays$RNA@data <- as.matrix(log2(tabula10X@assays$RNA@data + 1))
tabula10X <- FindVariableFeatures(tabula10X, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabula10X)
tabula10X <- ScaleData(tabula10X, features = all.genes)
tabula10X <- RunPCA(tabula10X, features = VariableFeatures(object = tabula10X))
tabula10X <- FindNeighbors(tabula10X, dims = 1:15)
tabula10X <- FindClusters(tabula10X, resolution = 0.4)
tabula10X <- RunUMAP(tabula10X, dims = 1:15, min.dist = 0.75)

FeaturePlot(tabula10X, features = c("logRNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (10X)") # 450 x 350

DimPlot(tabula10X, reduction = "umap", group.by = "tissue", pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (10X)") # height = 500


# Tabula 10X, absolute scaling + absolute transcripts projected
tabula10X@assays$RNA@data <- as.matrix(log2(tabula10X@assays$RNA@counts + 1))
tabula10X <- FindVariableFeatures(tabula10X, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabula10X)
tabula10X <- ScaleData(tabula10X, features = all.genes)
tabula10X <- RunPCA(tabula10X, features = VariableFeatures(object = tabula10X))
tabula10X <- FindNeighbors(tabula10X, dims = 1:15)
tabula10X <- FindClusters(tabula10X, resolution = 0.4)
tabula10X <- RunUMAP(tabula10X, dims = 1:15, min.dist = 0.75)

FeaturePlot(tabula10X, features = c("logRNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (10X)") # 450 x 350

DimPlot(tabula10X, reduction = "umap", group.by = "tissue", pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Tabula Muris (10X)") # height = 500