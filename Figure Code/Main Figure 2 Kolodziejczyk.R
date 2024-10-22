library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(data.table)
library(scran)
library(scRNAseq)
library(scater)
library(BBmisc)
library(cowplot)
library(tidyverse)

#Import Data
raw.data = as.sparse(read.csv(file = paste0("Kolodziejczyk2015/", "mesc.csv"), sep = ",", header = TRUE, row.names = 1))
celltype = as.data.frame(t(read.csv(file = paste0("Kolodziejczyk2015/", "celltype.csv"), sep = ",", header = TRUE, row.names = 1)))

mesc.seurat = CreateSeuratObject(counts = raw.data, project = "mesc")
mesc.seurat = AddMetaData(mesc.seurat, PercentageFeatureSet(mesc.seurat, pattern = "^ERCC-"), col.name = "ERCC")
mesc.seurat @meta.data$serum = celltype$ 'Serum-media'
mesc.seurat = subset(mesc.seurat, subset = ERCC < 80 & ERCC > 5 & nFeature_RNA > 1000 & nCount_RNA > 200)

ercc_spikecounts = mesc.seurat @meta.data$ERCC
filteredcounts = GetAssayData(object = mesc.seurat, assay = "RNA", slot = "counts")

#Absolute Scaling
scran.mesc = SingleCellExperiment(assays = list(counts = filteredcounts))
is.spike = grepl("^ERCC-", rownames(scran.mesc))
scran.mesc = splitAltExps(scran.mesc, ifelse(is.spike, "ERCC", "gene"))
scran.mesc = computeSpikeFactors(scran.mesc, "ERCC")
summary(sizeFactors(scran.mesc))
scran.mesc = logNormCounts(scran.mesc, log = FALSE)

mesc_ercc = CreateSeuratObject(counts = as.data.frame(x = assay(scran.mesc, "normcounts")))
mesc_ercc @meta.data$ERCC = ercc_spikecounts
mesc_ercc @meta.data$serum = mesc.seurat @meta.data$serum

mesc_ercc @assays$RNA @data = as.matrix(log2(mesc_ercc @assays$RNA @counts + 1))
mesc_ercc = FindVariableFeatures(mesc_ercc, selection.method = "vst", nfeatures = 10000)
all.genes = rownames(mesc_ercc)
mesc_ercc = ScaleData(mesc_ercc, features = all.genes)
mesc_ercc = RunPCA(mesc_ercc, features = VariableFeatures(object = mesc_ercc))

mesc_ercc = FindNeighbors(mesc_ercc, dims = 1: 10)
mesc_ercc = FindClusters(mesc_ercc, resolution = 0.2)
mesc_ercc = RunUMAP(mesc_ercc, dims = 1: 10)

Idents(mesc_ercc) = 'serum'
new.cluster.ids = c("2i", "Serum")
names(new.cluster.ids) = levels(mesc_ercc)
mesc_ercc = RenameIdents(mesc_ercc, new.cluster.ids)

saveRDS(mesc_ercc, "Kolodziejczyk2015/mesc_ercc")

#Global Scaling
scran.mesc = SingleCellExperiment(assays = list(counts = filteredcounts))
is.spike = grepl("^ERCC-", rownames(scran.mesc))
scran.mesc = splitAltExps(scran.mesc, ifelse(is.spike, "ERCC", "gene"))

normmat = as.data.frame(x = assay(scran.mesc, "counts"))
normmat = normmat[colnames(normmat) % in % colnames(mesc_ercc), ]
mesc_log = CreateSeuratObject(counts = normmat)
mesc_log = NormalizeData(object = mesc_log, normalization.method = "RC")
mesc_log = CreateSeuratObject(counts = as.data.frame(mesc_log @assays$RNA @data))
mesc_log = AddMetaData(mesc_log, PercentageFeatureSet(mesc_log, pattern = "^ERCC-"), col.name = "ERCC")
mesc_log @meta.data$serum = mesc.seurat @meta.data$serum

mesc_log @assays$RNA @data = as.matrix(log2(mesc_log @assays$RNA @counts + 1))
mesc_log = FindVariableFeatures(mesc_log, selection.method = "vst", nfeatures = 10000)
all.genes = rownames(mesc_log)
mesc_log = ScaleData(mesc_log, features = all.genes)
mesc_log = RunPCA(mesc_log, features = VariableFeatures(object = mesc_log))
mesc_log = FindNeighbors(mesc_log, dims = 1: 10)
mesc_log = FindClusters(mesc_log, resolution = 0.2)
mesc_log = RunUMAP(mesc_log, dims = 1: 10)

Idents(mesc_log) = 'serum'
new.cluster.ids = c("2i", "Serum")
names(new.cluster.ids) = levels(mesc_log)
mesc_log = RenameIdents(mesc_log, new.cluster.ids)

saveRDS(mesc_log, "Kolodziejczyk2015/mesc_log")

#Analysis
mesc_ercc = readRDS("Kolodziejczyk2015/mesc_ercc")
mesc_log = readRDS("Kolodziejczyk2015/mesc_log")
mesc_log @meta.data$nCount_logRNA = colSums(mesc_log @assays$RNA @data)
mesc_ercc @meta.data$celltype = Idents(mesc_ercc)
mesc_log @meta.data$celltype = Idents(mesc_log)

par(pty = "s")
dimcols = c("mediumaquamarine", "gold3")

DimPlot(mesc_ercc, reduction = "umap", pt.size = 1.5, cols = dimcols) +
  theme(aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
DimPlot(mesc_log, reduction = "umap", pt.size = 1.5, cols = dimcols) +
  theme(aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

FeaturePlot(mesc_ercc, features = c('nCount_RNA'), pt.size = 1.5) +
  theme(aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
        plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
FeaturePlot(mesc_log, features = c('nCount_logRNA'), pt.size = 1.5) +
  theme(aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
        plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

VlnPlot(mesc_ercc, features = c("nCount_RNA"), cols = dimcols, ncol = 1, pt.size = 0, y.max = (max(mesc_ercc$nCount_RNA) + 1000000)) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = "grey27", shape = 95) +
  scale_y_continuous(labels = function (x) format(x, scientific = TRUE)) +
  theme(aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
        plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0))

VlnPlot(mesc_log, features = c("nCount_logRNA"), cols = dimcols, ncol = 1, pt.size = 0, y.max = (max(mesc_log$nCount_RNA) + 1000000)) +
  ylab("log2 Transcript Count") + xlab("") + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = "grey27", shape = 95) +
  scale_y_continuous(labels = function (x) format(x, scientific = TRUE)) +
  theme(aspect.ratio = 1, text = element_text(size = 14), axis.text = element_text(size = 12),
        plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), axis.text.x = element_text(angle = 0))

#Transcriptome Curve
#Absolute scaling
ercc_cum = as.data.frame(t(as.data.frame(mesc_ercc @assays$RNA @data)))
ercc_cum = as.data.frame((scale(ercc_cum)))
ercc_cum = ercc_cum[, apply(ercc_cum, 2, function (x) !any(is.na(x)))]
ercc_cum[, "celltype"] = mesc_ercc @active.ident
ercc_cum = aggregate(.~celltype, ercc_cum, mean, na.action = na.omit)
rownames(ercc_cum) = ercc_cum$celltype
ercc_cum$celltype = NULL
ercc_cum = as.data.frame(t(ercc_cum))

ercc_tidy = ercc_cum % > %
  select(`2i`, Serum) % > %
  gather(key = "variable", value = "value")

ggplot(ercc_tidy, aes(x = value)) +
  stat_ecdf(aes(color = `variable`), geom = "step", size = 1.5) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_blank(),
        text = element_text(size = 14), axis.text = element_text(size = 12),
        plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  labs(title = "ERCC Normalized", y = "Cumulative Fraction", x = "Log2 Expression") +
  scale_color_manual(values = c("mediumaquamarine", "gold3")) +
  NoLegend()

#Global Scaling
log_cum = as.data.frame(t(as.data.frame(mesc_log @assays$RNA @data)))
log_cum = as.data.frame((scale(log_cum)))
log_cum = log_cum[, apply(log_cum, 2, function (x) !any(is.na(x)))]
log_cum[, "celltype"] = mesc_log @active.ident
log_cum = aggregate(.~celltype, log_cum, mean, na.action = na.omit)
rownames(log_cum) = log_cum$celltype
log_cum$celltype = NULL
log_cum = as.data.frame(t(log_cum))

log_tidy = log_cum % > %
  select(`2i`, Serum) % > %
  gather(key = "variable", value = "value")

ggplot(log_tidy, aes(x = value)) +
  stat_ecdf(aes(color = `variable`), geom = "step", size = 1.5) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_blank(),
        text = element_text(size = 14), axis.text = element_text(size = 12),
        plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  labs(title = "log Normalized", y = "Cumulative Fraction", x = "Log2 Expression") +
  scale_color_manual(values = c("mediumaquamarine", "gold3")) +
  NoLegend()

#Diff Expression Testing
levels(mesc_ercc)
mesc_ercc_markers = FindMarkers(mesc_ercc, ident = "Serum", slot = "data")
mesc_ercc_markers = mesc_ercc_markers[mesc_ercc_markers$p_val_adj < 0.05, ]
write.csv(mesc_ercc_markers, "Kolodziejczyk2015/mesc_ercc_diffgenes.csv")

levels(mesc_log)
mesc_log_markers = FindMarkers(mesc_log, ident = "Serum", slot = "data")
mesc_log_markers = mesc_log_markers[mesc_log_markers$p_val_adj < 0.05, ]
write.csv(mesc_log_markers, "Kolodziejczyk2015/mesc_log_diffgenes.csv")

#Gene Ontology Terms
ercc_go = read.csv("Kolodziejczyk2015/mesc_ercc_serum_Enrichr.csv")
ercc_go$logp = -log(ercc_go$Adjusted.P.value)
ercc_go$Name = sapply(X = strsplit(ercc_go$ï..Term, split = " \\("), FUN = "[", 1)
log_go = read.csv("Kolodziejczyk2015/mesc_log_serum_Enrichr.csv")
log_go$logp = -log(log_go$Adjusted.P.value)
log_go$Name = sapply(X = strsplit(log_go$ï..Term, split = " \\("), FUN = "[", 1)

ggplot(head(ercc_go, 10), aes(x = logp, y = reorder(Name, logp)), ) +
  geom_bar(stat = "identity", fill = "turquoise4") + theme_classic() +
  NoLegend() +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 0, unit = "pt"),
        text = element_text(size = 12), axis.text = element_text(size = 12),
        plot.title = element_blank(), aspect.ratio = 0.6,
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  xlim(0, 50) +
  labs(title = "", y = "", x = "-log p-Value")

ggplot(head(log_go, 10), aes(x = logp, y = reorder(Name, logp)), ) +
  geom_bar(stat = "identity", fill = "turquoise4") + theme_classic() +
  NoLegend() +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 0, unit = "pt"),
        text = element_text(size = 12), axis.text = element_text(size = 12),
        plot.title = element_blank(), aspect.ratio = 0.6,
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  xlim(0, 50) +
  labs(title = "", y = "", x = "-log p-Value")