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

tabula_skin <- readRDS("Tabula FACS Skin velocyto/skin.rds") # subset from Tabula10X/tabulafacsabs.rds
spliced <- readRDS("Tabula FACS Skin velocyto/spliced.rds")
spliced <- readRDS("Tabula FACS Skin velocyto/unspliced.rds")

cellnames <- colnames(tabula_skin)
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
scran.data <- splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))
altExpNames(scran.data)
scran.data <- computeSpikeFactors(scran.data, "ERCC")
summary(sizeFactors(scran.data))
scran.data <- logNormCounts(scran.data, log = FALSE)
spliced_ercc <- as.data.frame((x <- assay(scran.data, "normcounts")))
saveRDS(spliced_ercc, "Tabula FACS Skin velocyto/spliced_ercc.rds")

scran.data <- SingleCellExperiment(assays = list(counts = unspliced))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <- splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))
altExpNames(scran.data)
scran.data <- computeSpikeFactors(scran.data, "ERCC")
summary(sizeFactors(scran.data))
scran.data <- logNormCounts(scran.data, log = FALSE)
unspliced_ercc <- as.data.frame((x <- assay(scran.data, "normcounts")))
saveRDS(unspliced_ercc, "Tabula FACS Skin velocyto/unspliced_ercc.rds")

tabula_skin@meta.data$cellnames <- colnames(tabula_skin)
tabula_skin <- subset(tabula_skin, subset = cellnames %in% colnames(spliced_ercc))

tabula_skin[["spliced"]] <- CreateAssayObject(counts = spliced_ercc)
tabula_skin[["unspliced"]] <- CreateAssayObject(counts = unspliced_ercc)
rm(list = ls()[!ls() %in% c("tabula_skin")])

tabula_skin@assays$RNA@data <- as.matrix(log2(tabula_skin@assays$RNA@data + 1))
tabula_skin <- FindVariableFeatures(tabula_skin, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabula_skin)
tabula_skin <- ScaleData(tabula_skin, features = all.genes)
tabula_skin <- RunPCA(tabula_skin, features = VariableFeatures(object = tabula_skin))
tabula_skin <- FindNeighbors(tabula_skin, dims = 1:15)
tabula_skin <- FindClusters(tabula_skin, resolution = 0.4)
tabula_skin <- RunUMAP(tabula_skin, dims = 1:15, min.dist = 0.75)

# Panels G
FeaturePlot(tabula_skin, features = c("nCount_RNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Spliced Transcripts") # 500 x 300

FeaturePlot(tabula_skin, features = c("nCount_unspliced"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Unspliced Transcripts")

FeaturePlot(tabula_skin, features = c("Chd1"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Chd1")

FeaturePlot(tabula_skin, features = c("Myc"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Myc")

# Panels H
# Signatures
vircols <- viridis(10, direction = 1, option = "D")
FeaturePlot(tabula_skin, features = c("GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I")

FeaturePlot(tabula_skin, features = c("Hypertranscription"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Serum Hypertranscription")

FeaturePlot(tabula_skin, features = c("GO-RIBOSOME-BIOGENESIS"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-RIBOSOME-BIOGENESIS")

FeaturePlot(tabula_skin, features = c("GO-DNA-REPAIR"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-DNA-REPAIR")

# Panels A,D,F
DimPlot(tabula_skin, reduction = "umap", pt.size = 0.1, group.by = "celltype") +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

organcols <- rep("#00BF7D", 50)
library(ggridges)
figcelltype <- as.data.frame(tabula_skin@meta.data$celltype)
colnames(figcelltype) <- "CellType"
figcelltype[, "Transcripts"] <- tabula_skin@meta.data$logRNA

ggplot(figcelltype, aes(y = reorder(CellType, Transcripts, median), x = Transcripts)) +
  geom_density_ridges(scale = 2, fill = "#C77CFF") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme() +
  NoLegend() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(color = "black")) +
  labs(title = "", y = "", x = "log2 Transcripts")

Idents(tabula_skin) <- "celltype"
tabula_skin <- ReorderIdent(object = tabula_skin, var = "nCount_RNA", afxn = median)
DotPlot(tabula_skin, features = c(
  "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Top2b",
  "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
  "Rpl9", "Gapdh", "Actb", "Top2a", "Krt14", "Dkk3", "Fgf18"
)) +
  labs(title = "", y = "", x = "") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10.5)
  ) +
  NoLegend()

# Panels E
skin_curve <- as.data.frame(t(as.data.frame(tabula_skin@assays$RNA@data)))
skin_curve[, "celltype"] <- tabula_skin@meta.data$celltype
skin_curve <- aggregate(. ~ celltype, skin_curve, mean)
rownames(skin_curve) <- skin_curve$celltype
skin_curve$celltype <- NULL
skin_curve <- as.data.frame(t(skin_curve))
skin_curve[, "total"] <- rowSums(skin_curve)
skin_curve[, "genes"] <- rownames(skin_curve)
skin_curve <- skin_curve[skin_curve$total > 0, ]
skin_curve <- skin_curve[order(skin_curve$total), ]

library(tidyverse)
order.total <- order(skin_curve$total, decreasing = TRUE)
skin_curve$rank <- NA
skin_curve$rank[order.total] <- 1:nrow(skin_curve)

skin_tidy <- skin_curve %>%
  select(-c(total, genes)) %>%
  gather(key = "variable", value = "value", -rank)

skincols <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "#C77CFF")
ggplot(skin_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Stem Cell of Epidermis", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = skincols) +
  NoLegend()

skincols <- c("lightgrey", "lightgrey", "lightgrey", "#C77CFF", "lightgrey")
ggplot(skin_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Epidermal Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = skincols) +
  NoLegend()

skincols <- c("lightgrey", "#C77CFF", "lightgrey", "lightgrey", "lightgrey")
ggplot(skin_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Basal Cell of Epidermis", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = skincols) +
  NoLegend()

# Panels B-C
# Top2a Annotation
Top2a <- as.data.frame(tabula_skin@assays$RNA@counts["Top2a", ])
colnames(Top2a) <- "Top2a"
for (row in 1:nrow(Top2a)) {
  expr <- Top2a[row, "Top2a"]
  if (expr > 2.5) {
    Top2a[row, "status"] <- "Top2a+"
  } else {
    Top2a[row, "status"] <- "Top2a-"
  }
} # assign PGCs by Top2a expression
tabula_skin@meta.data$Top2a <- Top2a$status

DimPlot(tabula_skin, reduction = "umap", pt.size = 0.1, group.by = "Top2a", cols = c("lightgrey", "#C77CFF")) +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

VlnPlot(tabula_skin, features = c("nCount_RNA"), group.by = "Top2a", cols = c("lightgrey", "#C77CFF"), ncol = 1, pt.size = 0) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  labs(title = "")
