library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)
library(viridis)
library(monocle3)
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(parallelDist)
library(scran)
library(SingleCellExperiment)
library(data.table)
library(scRNAseq)
library(scater)

stomach <- readRDS("Busslinger2021/stomach.rds") # merged count matrices of pylorus/corpus
stomach@meta.data$noclusters <- 0
stomach <- AddMetaData(stomach, PercentageFeatureSet(stomach, pattern = "^ERCC."), col.name = "ERCC")
stomach <- subset(stomach, subset = nFeature_RNA > 1000 & nCount_RNA > 500 & ERCC > 0 & ERCC < 15)

filteredcounts <- GetAssayData(object = stomach, assay = "RNA", slot = "counts")
scran.stomach <- SingleCellExperiment(assays = list(counts = filteredcounts))
is.spike <- grepl("^ERCC.", rownames(scran.stomach))
scran.stomach <- splitAltExps(scran.stomach, ifelse(is.spike, "ERCC", "gene"))
scran.stomach <- computeSpikeFactors(scran.stomach, "ERCC")
summary(sizeFactors(scran.stomach))
scran.stomach <- logNormCounts(scran.stomach, log = FALSE)

stomach_ercc <- CreateSeuratObject(counts = as.data.frame(x = assay(scran.stomach, "normcounts")))
stomach_ercc@meta.data$noclusters <- 0
stomach_ercc <- subset(stomach_ercc, nCount_RNA < 250000)

stomach_ercc@assays$RNA@data <- as.matrix(log2(stomach_ercc@assays$RNA@counts + 1))
stomach_ercc <- FindVariableFeatures(stomach_ercc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(stomach_ercc)
stomach_ercc <- ScaleData(stomach_ercc, features = all.genes)
stomach_ercc <- RunPCA(stomach_ercc, features = VariableFeatures(object = stomach_ercc))
stomach_ercc <- FindNeighbors(stomach_ercc, dims = 1:15)
stomach_ercc <- FindClusters(stomach_ercc, resolution = 0.4)
stomach_ercc <- RunUMAP(stomach_ercc, dims = 1:15, min.dist = 0.75)

# Panels E
FeaturePlot(stomach_ercc, features = c("nCount_RNA"), pt.size = 1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Transcripts")

FeaturePlot(stomach_ercc, features = c("Mki67"), pt.size = 1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Mki67 (Proliferating Isthmus Cells)")

FeaturePlot(stomach_ercc, features = c("Muc6"), pt.size = 1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Muc6 (Mucous Gland Cell)")

FeaturePlot(stomach_ercc, features = c("Atp4b"), pt.size = 1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Atp4b (Parietal Cell)")

FeaturePlot(stomach_ercc, features = c("Anpep"), pt.size = 1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Anpep (Chief Cell)")

cluster.id <- c(
  "Mucous gland cell",
  "Surface mucous cell",
  "Surface mucous cell",
  "Mucous neck cell",
  "Neuroendocrine cell",
  "Surface mucous cell",
  "Chief cell",
  "Parietal cell",
  "Mast cell",
  "Isthmus cell",
  "Chief cell",
  "Tuft cell",
  "Immune cell"
)

names(cluster.id) <- levels(stomach_ercc)
stomach_ercc <- RenameIdents(stomach_ercc, cluster.id)
stomach_ercc@meta.data$celltype <- stomach_ercc@active.ident

# Panel A
DimPlot(stomach_ercc, reduction = "umap", pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 12),
    plot.title = element_blank(),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )
+NoLegend()

# Panel B
dimcols <- c(rep("darkorange3", 10))
VlnPlot(stomach_ercc, group.by = "celltype", features = c("nCount_RNA"), cols = dimcols, ncol = 1, pt.size = 0.1, sort = TRUE) +
  ylab("Transcripts") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.5, text = element_text(size = 14), axis.text = element_text(size = 12),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 45)
  ) # 600 x 350

# Panel C
stomach_ercc <- ReorderIdent(object = stomach_ercc, var = "nCount_RNA", afxn = mean)
DotPlot(stomach_ercc, features = c(
  "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Top2b",
  "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
  "Rpl9", "Gapdh", "Actb", "Top2a", "Muc6", "Mki67", "Anpep", "Atp4b"
)) +
  labs(title = "", y = "", x = "") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10.5)
  ) +
  NoLegend()

# Panels D
# Gene Curve
stomach_curve <- as.data.frame(t(as.data.frame(stomach_ercc@assays$RNA@data)))
stomach_curve[, "celltype"] <- stomach_ercc@meta.data$celltype
stomach_curve <- aggregate(. ~ celltype, stomach_curve, mean)
rownames(stomach_curve) <- stomach_curve$celltype
stomach_curve$celltype <- NULL
stomach_curve <- as.data.frame(t(stomach_curve))
stomach_curve[, "total"] <- rowSums(stomach_curve)
stomach_curve[, "genes"] <- rownames(stomach_curve)
stomach_curve <- stomach_curve[stomach_curve$total > 0, ]
stomach_curve <- stomach_curve[order(stomach_curve$total), ]

library(tidyverse)
order.total <- order(stomach_curve$`total`, decreasing = TRUE)
stomach_curve$rank <- NA
stomach_curve$rank[order.total] <- 1:nrow(stomach_curve)

stomach_tidy <- stomach_curve %>%
  dplyr::select(-c(total, genes)) %>%
  gather(key = "variable", value = "value", -rank)

stomachcols <- c(rep("lightgrey", 10))
ggplot(stomach_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Parietal Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = stomachcols) +
  geom_smooth(data = stomach_tidy[stomach_tidy$variable == "Parietal cell", ], aes(x = rank, y = value), color = "darkorange3", size = 1.2, alpha = 1, se = F) +
  NoLegend()

ggplot(stomach_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Isthmus Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = stomachcols) +
  geom_smooth(data = stomach_tidy[stomach_tidy$variable == "Isthmus cell", ], aes(x = rank, y = value), color = "darkorange3", size = 1.2, alpha = 1, se = F) +
  NoLegend()

ggplot(stomach_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Neuroendocrine Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = stomachcols) +
  geom_smooth(data = stomach_tidy[stomach_tidy$variable == "Neuroendocrine cell", ], aes(x = rank, y = value), color = "darkorange3", size = 1.2, alpha = 1, se = F) +
  NoLegend()

ggplot(stomach_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Mast Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = stomachcols) +
  geom_smooth(data = stomach_tidy[stomach_tidy$variable == "Mast cell", ], aes(x = rank, y = value), color = "darkorange3", size = 1.2, alpha = 1, se = F) +
  NoLegend()

ggplot(stomach_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Chief Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = stomachcols) +
  geom_smooth(data = stomach_tidy[stomach_tidy$variable == "Chief cell", ], aes(x = rank, y = value), color = "darkorange3", size = 1.2, alpha = 1, se = F) +
  NoLegend()
