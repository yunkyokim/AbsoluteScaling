library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(scran)
library(scRNAseq)
library(scater)
library(DoubletFinder)
library(BBmisc)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(VISION)
library(viridis)

# Data Processing
tabula_marrow <- readRDS("marrow.rds") # subset from Tabula FACS/tabulafacsabs.rds
spliced <- readRDS("Tabula Muris FACS BM velocyto/spliced.rds") # generated with velocyto pipeline
unspliced <- readRDS("Tabula Muris FACS BM velocyto/unspliced.rds")

# Preprocessing + ERCC-Normalization
cellnames <- colnames(tabula_marrow)
spliced <- spliced[rownames(spliced) %in% cellnames, ]
unspliced <- unspliced[rownames(unspliced) %in% cellnames, ]
spliced <- as.data.frame(t(spliced))
unspliced <- as.data.frame(t(unspliced))

rows <- rownames(spliced)
rows <- gsub("ERCC.", "ERCC-", rows)
rownames(spliced) <- rows
rows <- rownames(unspliced)
rows <- gsub("ERCC.", "ERCC-", rows)
rownames(unspliced) <- rows
unspliced[1:92, ] <- spliced[1:92, ]

scran.data <- SingleCellExperiment(assays = list(counts = spliced))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <-
  splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))
altExpNames(scran.data)
scran.data <- computeSpikeFactors(scran.data, "ERCC")
summary(sizeFactors(scran.data))
scran.data <- logNormCounts(scran.data, log = FALSE)
spliced_ercc <- as.data.frame((x <- assay(scran.data, "normcounts")))

scran.data <- SingleCellExperiment(assays = list(counts = unspliced))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <-
  splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))
altExpNames(scran.data)
scran.data <- computeSpikeFactors(scran.data, "ERCC")
summary(sizeFactors(scran.data))
scran.data <- logNormCounts(scran.data, log = FALSE)
unspliced_ercc <- as.data.frame((x <- assay(scran.data, "normcounts")))

tabula_marrow[["spliced"]] <-
  CreateAssayObject(counts = spliced_ercc)
tabula_marrow[["unspliced"]] <-
  CreateAssayObject(counts = unspliced_ercc)

tabula_marrow@assays$RNA@data <- as.matrix(log2(tabula_marrow@assays$RNA@data +
  1))
tabula_marrow <-
  FindVariableFeatures(tabula_marrow,
    selection.method = "vst",
    nfeatures = 2000
  )
all.genes <- rownames(tabula_marrow)
tabula_marrow <- ScaleData(tabula_marrow, features = all.genes)
tabula_marrow <-
  RunPCA(tabula_marrow, features = VariableFeatures(object = tabula_marrow))
tabula_marrow <- FindNeighbors(tabula_marrow, dims = 1:15)
tabula_marrow <- FindClusters(tabula_marrow, resolution = 0.4)
tabula_marrow <-
  RunUMAP(tabula_marrow, dims = 1:15, min.dist = 0.75)

# Figure Panels
# Panel A, 950 x 600
DimPlot(
  tabula_marrow,
  reduction = "umap",
  pt.size = 0.1,
  group.by = "celltype"
) +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

# Panel F, 500 x 300
FeaturePlot(tabula_marrow,
  features = c("nCount_RNA"),
  pt.size = 0.1
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Spliced Transcripts")

FeaturePlot(tabula_marrow,
  features = c("nCount_unspliced"),
  pt.size = 0.1
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Unspliced Transcripts")

FeaturePlot(tabula_marrow,
  features = c("Chd1"),
  pt.size = 0.1
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Chd1")

# Panel E, 500 x 300
vircols <- viridis(10, direction = 1, option = "D")
FeaturePlot(
  tabula_marrow,
  features = c("GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I"),
  pt.size = 0.1,
  cols = vircols
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I")

FeaturePlot(
  tabula_marrow,
  features = c("Hypertranscription"),
  pt.size = 0.1,
  cols = vircols
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Serum Hypertranscription")

FeaturePlot(
  tabula_marrow,
  features = c("GO-RIBOSOME-BIOGENESIS"),
  pt.size = 0.1,
  cols = vircols
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "GO-RIBOSOME-BIOGENESIS")


library(ggridges)
figcelltype <- as.data.frame(tabula_marrow@meta.data$celltype)
colnames(figcelltype) <- "CellType"
figcelltype[, "Transcripts"] <- tabula_marrow@meta.data$logRNA

# Panel B, 650 x 1100
ggplot(figcelltype, aes(
  y = reorder(CellType, Transcripts, median), x =
    Transcripts
)) +
  geom_density_ridges(scale = 2, fill = "#35A2FF") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme() +
  NoLegend() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_text(color = "black")
  ) +
  labs(title = "", y = "", x = "log2 Transcripts") # 450 x 500

Idents(tabula_marrow) <- "celltype"
tabula_marrow <- ReorderIdent(object = tabula_marrow, var = "nCount_RNA", afxn = median)
# Panel D, 880 x 700
DotPlot(
  tabula_marrow,
  features = c(
    "Chd1",
    "Hira",
    "Ino80",
    "Kat5",
    "Myc",
    "Top2b",
    "Atm",
    "Parp1",
    "Rps2",
    "Rpl3",
    "Rpl4",
    "Rpl6",
    "Rpl9",
    "Gapdh",
    "Actb",
    "Cd34",
    "Cd19",
    "Cd74",
    "Cd4"
  )
) +
  labs(title = "", y = "", x = "") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(size = 10.5)
  ) +
  NoLegend()

# Gene Curves
marrow_curve <- as.data.frame(t(as.data.frame(tabula_marrow@assays$RNA@data)))
marrow_curve <- as.data.frame(t(as.data.frame(tabula_marrow@assays$RNA@data)))
marrow_curve[, "celltype"] <- tabula_marrow@meta.data$celltype
marrow_curve <- aggregate(. ~ celltype, marrow_curve, mean)
rownames(marrow_curve) <- marrow_curve$celltype
marrow_curve$celltype <- NULL
marrow_curve <- as.data.frame(t(marrow_curve))
marrow_curve[, "total"] <- rowSums(marrow_curve)
marrow_curve[, "genes"] <- rownames(marrow_curve)
marrow_curve <- marrow_curve[marrow_curve$total > 0, ]
marrow_curve <- marrow_curve[order(marrow_curve$total), ]

library(tidyverse)
order.total <- order(marrow_curve$total, decreasing = TRUE)
marrow_curve$rank <- NA
marrow_curve$rank[order.total] <- 1:nrow(marrow_curve)

marrow_tidy <- marrow_curve %>%
  select(-c(total, genes)) %>%
  gather(key = "variable", value = "value", -rank)

# Panel C Left, 600 x 350
marrowcols <- c(rep("lightgrey", 21), "#35A2FF")
ggplot(marrow_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 7,
      r = 0,
      b = 7,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Slamf1-positive multipotent progenitor cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = marrowcols) +
  NoLegend()

# Panel C Middle, 600 x 350
marrowcols <- c(rep("lightgrey", 15), "#35A2FF", rep("lightgrey", 6))
ggplot(marrow_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 7,
      r = 0,
      b = 7,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Monocyte", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = marrowcols) +
  NoLegend()

# Panel C Right, 600 x 350
marrowcols <- c(rep("lightgrey", 16), "#35A2FF", rep("lightgrey", 5))
ggplot(marrow_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 7,
      r = 0,
      b = 7,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Naive B Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = marrowcols) +
  NoLegend()
