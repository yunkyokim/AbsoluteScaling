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

tabula_intestine <- readRDS("Tabula FACS Intestine velocyto/intestine.rds") # subset from Tabula10X/tabulafacsabs.rds
spliced <- readRDS("Tabula FACS Intestine velocyto/spliced.rds")
spliced <- readRDS("Tabula FACS Intestine velocyto/unspliced.rds")

cellnames <- colnames(tabula_intestine)
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
saveRDS(spliced_ercc, "Tabula FACS Intestine velocyto/spliced_ercc.rds")

scran.data <- SingleCellExperiment(assays = list(counts = unspliced))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <- splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))
altExpNames(scran.data)
scran.data <- computeSpikeFactors(scran.data, "ERCC")
summary(sizeFactors(scran.data))
scran.data <- logNormCounts(scran.data, log = FALSE)
unspliced_ercc <- as.data.frame((x <- assay(scran.data, "normcounts")))
saveRDS(unspliced_ercc, "Tabula FACS Intestine velocyto/unspliced_ercc.rds")

tabula_intestine@meta.data$cellnames <- colnames(tabula_intestine)
tabula_intestine <- subset(tabula_intestine, subset = cellnames %in% colnames(spliced_ercc))

tabula_intestine[["spliced"]] <- CreateAssayObject(counts = spliced_ercc)
tabula_intestine[["unspliced"]] <- CreateAssayObject(counts = unspliced_ercc)
rm(list = ls()[!ls() %in% c("tabula_intestine")])

tabula_intestine@assays$RNA@data <- as.matrix(log2(tabula_intestine@assays$RNA@data + 1))
tabula_intestine <- FindVariableFeatures(tabula_intestine, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabula_intestine)
tabula_intestine <- ScaleData(tabula_intestine, features = all.genes)
tabula_intestine <- RunPCA(tabula_intestine, features = VariableFeatures(object = tabula_intestine))
tabula_intestine <- FindNeighbors(tabula_intestine, dims = 1:15)
tabula_intestine <- FindClusters(tabula_intestine, resolution = 0.4)
tabula_intestine <- RunUMAP(tabula_intestine, dims = 1:15, min.dist = 0.75)

# Panels G
FeaturePlot(tabula_intestine, features = c("nCount_RNA"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Spliced Transcripts") # 500 x 300

FeaturePlot(tabula_intestine, features = c("nCount_unspliced"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Unspliced Transcripts")

FeaturePlot(tabula_intestine, features = c("Chd1"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Chd1")

FeaturePlot(tabula_intestine, features = c("Myc"), pt.size = 0.1) +
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
FeaturePlot(tabula_intestine, features = c("GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I") # 500 x 300

FeaturePlot(tabula_intestine, features = c("Hypertranscription"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Serum Hypertranscription") # 500 x 300

FeaturePlot(tabula_intestine, features = c("GO-RIBOSOME-BIOGENESIS"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-RIBOSOME-BIOGENESIS") # 500 x 300

FeaturePlot(tabula_intestine, features = c("GO-DNA-REPAIR"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-DNA-REPAIR") # 500 x 300

# Panels A,D,F
DimPlot(tabula_intestine, reduction = "umap", pt.size = 0.1, group.by = "celltype") +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

organcols <- rep("#00BF7D", 50)
library(ggridges)
figcelltype <- as.data.frame(tabula_intestine@meta.data$celltype)
colnames(figcelltype) <- "CellType"
figcelltype[, "Transcripts"] <- tabula_intestine@meta.data$logRNA

ggplot(figcelltype, aes(y = reorder(CellType, Transcripts, median), x = Transcripts)) +
  geom_density_ridges(scale = 2, fill = "#00BF7D") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme() +
  NoLegend() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(color = "black")) +
  labs(title = "", y = "", x = "log2 Transcripts") # 500 x 250

Idents(tabula_intestine) <- "celltype"
tabula_intestine <- ReorderIdent(object = tabula_intestine, var = "nCount_RNA", afxn = median)
DotPlot(tabula_intestine, features = c(
  "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Top2b",
  "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
  "Rpl9", "Gapdh", "Actb", "Lgr5", "Krt20", "Atoh1", "Chga"
)) +
  labs(title = "", y = "", x = "") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10.5)
  ) +
  NoLegend()

# Panels E
intestine_curve <- as.data.frame(t(as.data.frame(tabula_intestine@assays$RNA@data)))
intestine_curve[, "celltype"] <- tabula_intestine@meta.data$celltype
intestine_curve <- aggregate(. ~ celltype, intestine_curve, mean)
rownames(intestine_curve) <- intestine_curve$celltype
intestine_curve$celltype <- NULL
intestine_curve <- as.data.frame(t(intestine_curve))
intestine_curve[, "total"] <- rowSums(intestine_curve)
intestine_curve[, "genes"] <- rownames(intestine_curve)
intestine_curve <- intestine_curve[intestine_curve$total > 0, ]
intestine_curve <- intestine_curve[order(intestine_curve$total), ]

library(tidyverse)
order.total <- order(intestine_curve$total, decreasing = TRUE)
intestine_curve$rank <- NA
intestine_curve$rank[order.total] <- 1:nrow(intestine_curve)

intestine_tidy <- intestine_curve %>%
  select(-c(total, genes)) %>%
  gather(key = "variable", value = "value", -rank)

intestinecols <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "#00BF7D")
ggplot(intestine_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Goblet Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = intestinecols) +
  NoLegend()

intestinecols <- c("lightgrey", "lightgrey", "lightgrey", "#00BF7D", "lightgrey")
ggplot(intestine_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Epithelial Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = intestinecols) +
  NoLegend()

intestinecols <- c("#00BF7D", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
ggplot(intestine_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Enteroendocrine Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = intestinecols) +
  NoLegend()

# Panels B-C
# Lgr5 Annotation
Lgr5 <- as.data.frame(tabula_intestine@assays$RNA@counts["Lgr5", ])
colnames(Lgr5) <- "Lgr5"
for (row in 1:nrow(Lgr5)) {
  expr <- Lgr5[row, "Lgr5"]
  if (expr > 2.5) {
    Lgr5[row, "status"] <- "Lgr5+"
  } else {
    Lgr5[row, "status"] <- "Lgr5-"
  }
} # assign by Lgr5 expression
tabula_intestine@meta.data$lgr5 <- Lgr5$status

DimPlot(tabula_intestine, reduction = "umap", pt.size = 0.1, group.by = "lgr5", cols = c("lightgrey", "#00BF7D")) +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

VlnPlot(tabula_intestine, features = c("nCount_RNA"), group.by = "lgr5", cols = c("lightgrey", "#00BF7D"), ncol = 1, pt.size = 0) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
