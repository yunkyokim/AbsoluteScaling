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


tabulafacs <- readRDS("Tabula FACS/tabulafacs_vision.rds")
tabulafacs@meta.data$logRNA <- log2(tabulafacs@meta.data$nCount_RNA + 1)
cols <- read.csv(file = "Tabula FACS/tabulafacscols.csv")

SIGcomparison <- as.data.frame(t(tabulafacs@assays$vision@data))
sigscores <- as.data.frame(t(tabulafacs@assays$hyper@data))
SIGcomparison[, "Hypertranscription"] <- sigscores$Hypertranscription
SIGcomparison[, "RNA_Counts"] <- tabulafacs@meta.data$logRNA

sig <- colnames(SIGcomparison)
sig <- gsub(
  x = sig,
  pattern = "\\-",
  replacement = "_"
)
colnames(SIGcomparison) <- sig
rm(sig)

library(MASS)
library(ggplot2)
library(viridis)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

ggplot(SIGcomparison, aes(x = RNA_Counts, y = GO_DNA_REPAIR, colour = get_density(RNA_Counts, GO_DNA_REPAIR, n = 200))) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm", se = TRUE) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  NoLegend() +
  labs(y = "GO_DNA_REPAIR", x = "Log2 Transcripts") # 600 x 500

ggplot(SIGcomparison, aes(x = RNA_Counts, y = GO_CHROMATIN_ORGANIZATION, colour = get_density(RNA_Counts, GO_CHROMATIN_ORGANIZATION, n = 200))) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm", se = TRUE) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  NoLegend() +
  labs(y = "GO_CHROMATIN_ORGANIZATION", x = "Log2 Transcripts") # 600 x 500

ggplot(SIGcomparison, aes(x = RNA_Counts, y = GO_CELL_DIVISION, colour = get_density(RNA_Counts, GO_CELL_DIVISION, n = 200))) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm", se = TRUE) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  NoLegend() +
  labs(y = "GO_CELL_DIVISION", x = "Log2 Transcripts") # 600 x 500

# Gene Comparison, organ level
GENEcomparison <- readRDS("Tabula FACS/tabulafacs_genecomparison.rds")

GENEcomparison[, "organcell"] <- rownames(GENEcomparison)
GENEcomparison[, "Organ"] <- sapply(X = strsplit(GENEcomparison$organcell, split = "/"), FUN = "[", 2)
GENEcomparison[, "Celltype"] <- sapply(X = strsplit(GENEcomparison$organcell, split = "/"), FUN = "[", 1)

GENEcomparison[, "colors"] <- cols$colors[match(GENEcomparison$Organ, cols$Organs)]
GENEcomparison$organcell <- paste(GENEcomparison$Celltype, " (", GENEcomparison$Organ, ")", sep = "")
rownames(GENEcomparison) <- GENEcomparison$organcell

ggplot(GENEcomparison, aes(x = RNA_Counts, y = Myc)) +
  geom_point(aes(color = colors), size = 2) +
  geom_text_repel(aes(label = ifelse(RNA_Counts > 5 & Myc > 3, as.character(Celltype), ""), color = colors),
    size = 3, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(y = "Myc", x = "Log2 Transcripts") +
  NoLegend()

ggplot(GENEcomparison, aes(x = RNA_Counts, y = Yap1)) +
  geom_point(aes(color = colors), size = 2) +
  geom_text_repel(aes(label = ifelse(RNA_Counts > 3 & Yap1 > 2, as.character(Celltype), ""), color = colors),
    size = 3, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(y = "Yap", x = "Log2 Transcripts") +
  NoLegend()

ggplot(GENEcomparison, aes(x = RNA_Counts, y = Mtor)) +
  geom_point(aes(color = colors), size = 2) +
  geom_text_repel(aes(label = ifelse(RNA_Counts > 6 & Mtor > 6, as.character(Celltype), ""), color = colors),
    size = 4, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(y = "Mtor", x = "Log2 Transcripts") +
  NoLegend()

ggplot(GENEcomparison, aes(x = RNA_Counts, y = Ets1)) +
  geom_point(aes(color = colors), size = 2) +
  geom_text_repel(aes(label = ifelse(RNA_Counts < 5 & Ets1 > 5, as.character(Celltype), ""), color = colors),
    size = 3, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(y = "Ets1", x = "Log2 Transcripts") +
  NoLegend()

ggplot(GENEcomparison, aes(x = RNA_Counts, y = Grap)) +
  geom_point(aes(color = colors), size = 2) +
  geom_text_repel(aes(label = ifelse(RNA_Counts < 5 & Grap > 5, as.character(Celltype), ""), color = colors),
    size = 3, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(y = "Grap", x = "Log2 Transcripts") +
  NoLegend()

ggplot(GENEcomparison, aes(x = RNA_Counts, y = Gbp9)) +
  geom_point(aes(color = colors), size = 2) +
  geom_text_repel(aes(label = ifelse(RNA_Counts < 5 & Gbp9 > 5, as.character(Celltype), ""), color = colors),
    size = 3, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(y = "Gbp9", x = "Log2 Transcripts") +
  NoLegend()

# CHEA
enrichrCHEA <- read.csv("Tabula FACS/enrichrCHEA.csv", fill = TRUE) # performed externally using highly correlated genes from GENEtest
enrichrCHEA$logp <- -log(enrichrCHEA$Adjusted.P.value)
enrichrCHEA$Name <- sapply(X = strsplit(enrichrCHEA$Term, split = " \\("), FUN = "[", 1)

ggplot(head(enrichrCHEA, 20), aes(x = logp, y = reorder(Name, logp)), ) +
  geom_bar(stat = "identity", fill = "maroon4") +
  theme_classic() +
  NoLegend() +
  theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0, unit = "pt"),
    text = element_text(size = 12), axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5), aspect.ratio = 1,
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  xlim(0, 230) +
  labs(title = "ChEA Transcription Factors", y = "", x = "-log p-Value") +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 45)) # 800 x 600
