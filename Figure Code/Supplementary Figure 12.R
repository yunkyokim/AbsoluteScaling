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

# Independently generated adult BM dataset, see methods
multi <- readRDS("Kim et al 10X BM 2022/bm_seurat_vision.rds")

# Panels
DimPlot(multi, reduction = "umap", pt.size = 1, label = T) +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

multi@assays$RNA@data <- as.matrix(log2(multi@assays$RNA@counts + 1))
multi <- ReorderIdent(object = multi, var = "nCount_RNA", afxn = median)
DotPlot(multi, features = c(
  "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Top2b",
  "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
  "Rpl9", "Gapdh", "Actb", "Cd34", "Cd19", "Cd74", "Cd4"
), cols = "PRGn") +
  labs(title = "", y = "", x = "") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10.5)
  ) +
  NoLegend()

multi@meta.data$celltype <- Idents(multi)
marrow_curve <- as.data.frame(t(as.data.frame(multi@assays$RNA@data)))
marrow_curve[, "celltype"] <- multi@meta.data$celltype
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

VlnPlot(multi, c("nCount_RNA", "nFeature_RNA"), group.by = "celltype")

marrowcols <- c(rep("lightgrey", 12), "#35A2FF")
ggplot(marrow_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Multipotent progenitor cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = marrowcols) +
  NoLegend()

marrowcols <- c("#35A2FF", rep("lightgrey", 12))
ggplot(marrow_tidy, aes(x = rank, y = value)) +
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
  labs(title = "Erythroblast", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = marrowcols) +
  NoLegend()

library(ggridges)
figcelltype <- as.data.frame(multi@meta.data$celltype)
colnames(figcelltype) <- "CellType"
multi@meta.data$log_RNA <- log(multi@meta.data$nCount_RNA + 1)
figcelltype[, "Transcripts"] <- multi@meta.data$log_RNA

ggplot(figcelltype, aes(y = reorder(CellType, Transcripts, median), x = Transcripts)) +
  geom_density_ridges(scale = 2, fill = "#35A2FF") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme() +
  NoLegend() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(color = "black")) +
  labs(title = "", y = "", x = "log2 Transcripts")

FeaturePlot(multi, features = c("log_RNA"), pt.size = 1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Log Transcripts")

FeaturePlot(multi, features = c("Chd1"), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Chd1")

# mSigDB GO BP Signatures
multi@assays$RNA@data <- as.matrix(log2(multi@assays$RNA@data + 1))
signature.list <- ("Tabula FACS/c5.bp.v7.1.symbols.gmt")
multimat <- as.sparse(multi@assays$RNA@counts)
multi.vision <- Vision(multimat, signatures = signature.list, projection_methods = NULL, pool = FALSE, sig_gene_threshold = 0.05, sig_norm_method = "none")
multi.vision <- calcSignatureScores(multi.vision, sig_norm_method = "none", sig_gene_importance = FALSE)

sigscores <- t(as.data.frame(getSignatureScores(multi.vision)))
multi[["vision"]] <- CreateAssayObject(counts = sigscores)

serum_diff <- read.csv("Tabula FACS/mesc_serum_sig.csv", header = FALSE)
serum_diff[1, 1] <- "TagIn"
serum_list <- (serum_diff$V2)
names(serum_list) <- serum_diff$V1

hypertranscription <- createGeneSignature("Hypertranscription", sigData = serum_list)
multi.vision <- Vision(multimat, signatures = c(hypertranscription), projection_methods = NULL, pool = FALSE, sig_norm_method = "none")
multi.vision <- calcSignatureScores(multi.vision, sig_norm_method = "znorm_rows", sig_gene_importance = FALSE)

sigscores <- t(as.data.frame(getSignatureScores(multi.vision)))
multi[["hyper"]] <- CreateAssayObject(counts = sigscores)

vircols <- viridis(10, direction = 1, option = "D")
FeaturePlot(multi, features = c("Hypertranscription"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "Serum Hypertranscription") # 500 x 300

FeaturePlot(multi, features = c("GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I") # 500 x 300

FeaturePlot(multi, features = c("GO-RIBOSOME-BIOGENESIS"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-RIBOSOME-BIOGENESIS") # 500 x 300

FeaturePlot(multi, features = c("GO-CHROMATIN-REMODELING"), pt.size = 0.1, cols = vircols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  ) +
  labs(title = "GO-CHROMATIN-REMODELING") # 500 x 300
