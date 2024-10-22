library(Seurat)
library(ggplot2)
library(BBmisc)
library(data.table)
library(dplyr)
library(scran)
library(scater)
library(tidyverse)

# Panel B - C
# Cell Hashing
hashing <- readRDS("Cell Hashing/hashing.subset")
hashing <- subset(hashing, subset = nFeature_RNA > 150 & nCount_RNA > 100 & nCount_RNA < 3000)

hashing @assays$RNA @data <- as.matrix(log2(hashing @assays$RNA @counts + 1))
hashingfc <- FindMarkers(hashing, "Doublet", "Singlet")

ribo_house <- c("RPL5", "RPL8", "RPL9", "RPL10A", "RPL11", "RPL14", "RPS5", "RPS6", "RPL13", "ACTB", "B2M", "GAPDH", "ALDOA", "PPIA")
DotPlot(hashing, features = ribo_house) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.35,
    axis.text = element_text(size = 13),
    plot.title = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

dotplot_fc <- data.frame()
for (gene in ribo_house) {
  averageexp <- as.data.frame(AverageExpression(hashing, features = gene))
  dotplot_fc[gene, "fc"] <- (log2(averageexp[1, 2] / averageexp[1, 1]))
}

# Mean RPL / RPS Expression
genelist <- rownames(hashing @assays$RNA @data)
ribolist <- c(genelist[grep(c("^(?=.*RPL)(?!.*MRPL)"), genelist, perl = TRUE)], genelist[grep(c("^(?=.*RPS)(?!.*MRPS)"), genelist, perl = TRUE)])
meanribo <- colMeans(hashing @assays$RNA @data[ribolist, ], na.rm = TRUE)
hashing[["ribo"]] <- meanribo

vlncol <- c("deepskyblue4", "deepskyblue3")
VlnPlot(hashing, features = c("ribo"), group.by = "orig.ident", cols = vlncol, ncol = 1, pt.size = 0) +
  ylab("log2 Expression") +
  xlab("") +
  NoLegend() +
  ggtitle("") +
  stat_summary(fun.y = median, geom = "point", size = 15, colour = "grey27", shape = 95) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# Demuxlet
demuxlet <- readRDS("Demuxlet/demuxlet.subset")
demuxlet @assays$RNA @data <- as.matrix(log2(demuxlet @assays$RNA @counts + 1))
Idents(demuxlet) <- "status"
demuxletfc <- FindMarkers(demuxlet, "Doublet", "Singlet")

ribo_house <- c("RPL5", "RPL8", "RPL9", "RPL10A", "RPL11", "RPL14", "RPS5", "RPS6", "RPL13", "ACTB", "B2M", "GAPDH", "ALDOA", "PPIA")
DotPlot(demuxlet, features = ribo_house) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.35,
    axis.text = element_text(size = 13),
    plot.title = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

dotplot_fc <- data.frame()
for (gene in ribo_house) {
  averageexp <- as.data.frame(AverageExpression(demuxlet, features = gene))
  dotplot_fc[gene, "fc"] <- ((averageexp[1, 2] / averageexp[1, 1]))
}

# Mean RPL / RPS Expression
genelist <- rownames(demuxlet @assays$RNA @data)
ribolist <- c(genelist[grep(c("^(?=.*RPL)(?!.*MRPL)"), genelist, perl = TRUE)], genelist[grep(c("^(?=.*RPS)(?!.*MRPS)"), genelist, perl = TRUE)])
meanribo <- colMeans(demuxlet @assays$RNA @data[ribolist, ], na.rm = TRUE)
demuxlet[["ribo"]] <- meanribo

vlncol <- c("deepskyblue4", "deepskyblue3")
VlnPlot(demuxlet, features = c("ribo"), group.by = "status", cols = vlncol, ncol = 1, pt.size = 0) +
  ylab("log2 Expression") +
  xlab("") +
  NoLegend() +
  ggtitle("") +
  stat_summary(fun.y = median, geom = "point", size = 15, colour = "grey27", shape = 95) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) # 400 x 300
