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

#setwd("current_dir/")

#Data import, files available upon request
#Import normalized Seurat objects
mesc_ercc = readRDS("mesc_ercc.rds")
mesc_log = readRDS("mesc_log.rds")

mesc_log@meta.data$nCount_logRNA = colSums(mesc_log@assays$RNA@data)
mesc_ercc@meta.data$celltype = Idents(mesc_ercc)
mesc_log@meta.data$celltype = Idents(mesc_log)


#Figure Panels
par(pty = "s")
dimcols = c("mediumaquamarine", "gold3")
#Panel A, 600 x 350
DimPlot(mesc_log,
        reduction = "umap",
        pt.size = 1.5,
        cols = dimcols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  )

#Panel B, 600 x 350
DimPlot(mesc_ercc,
        reduction = "umap",
        pt.size = 1.5,
        cols = dimcols) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  )

#Panel C, 600 x 350
FeaturePlot(mesc_log,
            features = c('nCount_logRNA'),
            pt.size = 1.5) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  )

#Panel D, 600 x 350
FeaturePlot(mesc_ercc,
            features = c('nCount_RNA'),
            pt.size = 1.5) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  )

#Panel E, 600 x 350
VlnPlot(
  mesc_log,
  features = c("nCount_logRNA"),
  cols = dimcols,
  ncol = 1,
  pt.size = 0,
  y.max = (max(mesc_log$nCount_RNA) + 1000000)
) +
  ylab("log2 Transcript Count") + xlab("") + NoLegend() +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 15,
    colour = "grey27",
    shape = 95
  ) +
  #stat_compare_means(comparisons = list(c("2i", "Serum")), label = "p.signif")+
  scale_y_continuous(
    labels = function(x)
      format(x, scientific = TRUE)
  ) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0)
  ) # 600 x 350

#Calulate fold-change, LOG
mesc_log_transcripts = data.frame(
  transcripts = mesc_log@meta.data$nCount_logRNA,
  identity = mesc_log@meta.data$celltype
)
mesc_log_transcripts = aggregate(. ~ identity, mesc_log_transcripts, median)
rownames(mesc_log_transcripts) = mesc_log_transcripts$identity
mesc_log_transcripts["Serum", "transcripts"] / mesc_log_transcripts["2i", "transcripts"]

#Panel D, 600 x 350
VlnPlot(
  mesc_ercc,
  features = c("nCount_RNA"),
  cols = dimcols,
  ncol = 1,
  pt.size = 0,
  y.max = (max(mesc_ercc$nCount_RNA) + 1000000)
) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 15,
    colour = "grey27",
    shape = 95
  ) +
  scale_y_continuous(
    labels = function(x)
      format(x, scientific = TRUE)
  ) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0)
  )


mesc_ercc_transcripts = data.frame(
  transcripts = mesc_ercc@meta.data$nCount_RNA,
  identity = mesc_ercc@meta.data$celltype
)
mesc_ercc_transcripts = aggregate(. ~ identity, mesc_ercc_transcripts, median)
rownames(mesc_ercc_transcripts) = mesc_ercc_transcripts$identity
mesc_ercc_transcripts["Serum", "transcripts"] / mesc_ercc_transcripts["2i", "transcripts"]


#Distribution Function
#Global Scaling
log_cum = as.data.frame(t(as.data.frame(mesc_log@assays$RNA@data)))
log_cum = as.data.frame((scale(log_cum)))
log_cum = log_cum[, apply(log_cum, 2, function(x)
  ! any(is.na(x)))]
log_cum[, "celltype"] = mesc_log@active.ident
log_cum = aggregate(. ~ celltype, log_cum, mean, na.action = na.omit)
rownames(log_cum) = log_cum$celltype
log_cum$celltype = NULL
log_cum = as.data.frame(t(log_cum))

log_tidy = log_cum %>%
  select(`2i`, Serum) %>%
  gather(key = "variable", value = "value")

#Panel H, 600 x 350
ggplot(log_tidy, aes(x = value)) +
  stat_ecdf(aes(color = `variable`), geom = "step", size = 1.5) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
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
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "log Normalized", y = "Cumulative Fraction", x = "Log2 Expression") +
  scale_color_manual(values = c("mediumaquamarine", "gold3")) +
  NoLegend()

#Absolute Scaling
ercc_cum = as.data.frame(t(as.data.frame(mesc_ercc@assays$RNA@data)))
ercc_cum = as.data.frame((scale(ercc_cum)))
ercc_cum = ercc_cum[, apply(ercc_cum, 2, function(x)
  ! any(is.na(x)))]
ercc_cum[, "celltype"] = mesc_ercc@active.ident
ercc_cum = aggregate(. ~ celltype, ercc_cum, mean, na.action = na.omit)
rownames(ercc_cum) = ercc_cum$celltype
ercc_cum$celltype = NULL
ercc_cum = as.data.frame(t(ercc_cum))

ercc_tidy = ercc_cum %>%
  select(`2i`, Serum) %>%
  gather(key = "variable", value = "value")

#Panel H, 600 x 350
ggplot(ercc_tidy, aes(x = value)) +
  stat_ecdf(aes(color = `variable`), geom = "step", size = 1.5) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
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
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "ERCC Normalized", y = "Cumulative Fraction", x = "Log2 Expression") +
  scale_color_manual(values = c("mediumaquamarine", "gold3")) +
  NoLegend()

#Differential Expression Testing
mesc_log_markers = FindMarkers(mesc_log, ident.1 = "Serum", slot = "data")
mesc_log_markers = mesc_log_markers[mesc_log_markers$p_val_adj < 0.05,]
write.csv(mesc_log_markers, "mesc_log_diffgenes.csv")

mesc_ercc_markers = FindMarkers(mesc_ercc, ident.1 = "Serum", slot = "data")
mesc_ercc_markers = mesc_ercc_markers[mesc_ercc_markers$p_val_adj < 0.05,]
write.csv(mesc_ercc_markers, "mesc_ercc_diffgenes.csv")

#Gene Ontology Terms
log_go = read.csv("mesc_log_serum_Enrichr.csv")
log_go$logp = -log(log_go$Adjusted.P.value)
log_go$Name = sapply(X = strsplit(log_go$ï..Term, split = " \\("), FUN = "[", 1)

ercc_go = read.csv("mesc_ercc_serum_Enrichr.csv")
ercc_go$logp = -log(ercc_go$Adjusted.P.value)
ercc_go$Name = sapply(X = strsplit(ercc_go$ï..Term, split = " \\("), FUN = "[", 1)

#Panel I, 800 x 250
ggplot(head(log_go, 10), aes(x = logp, y = reorder(Name, logp)),) +
  geom_bar(stat = "identity",  fill = "turquoise4") + theme_classic() +
  NoLegend() +
  theme(
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 0,
      unit = "pt"
    ),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_blank(),
    aspect.ratio = 0.6,
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  xlim(0, 50) +
  labs(title = "", y = "", x = "-log p-Value")

#Panel J, 800 x 250
ggplot(head(ercc_go, 10), aes(x = logp, y = reorder(Name, logp)),) +
  geom_bar(stat = "identity",  fill = "turquoise4") + theme_classic() +
  NoLegend() +
  theme(
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 0,
      unit = "pt"
    ),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_blank(),
    aspect.ratio = 0.6,
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  xlim(0, 50) +
  labs(title = "", y = "", x = "-log p-Value") #800 x 250
