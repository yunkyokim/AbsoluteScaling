library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(scran)
library(scater)
library(BBmisc)
library(tidyverse)

#setwd("current_dir/")
#Data import, SORT-seq or CEL-seq datasets
scattercol = c("seagreen4", "seagreen3", "seagreen2", "seagreen1")
umi_raw = read.csv(
  "UMI Gene Count.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
)
ercc_raw = read.csv(
  "ERCC Gene Count.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
)
barcodes = read.csv("barcodes.csv",
                    sep = ",",
                    header = TRUE,
                    row.names = 1)


#UMI Global-scaling by relative counts
UMI_Glo2 <- CreateSeuratObject(counts = umi_raw)
UMI_Glo2 = AddMetaData(UMI_Glo2,
                       metadata = barcodes$mRNA_amount,
                       col.name = "mRNAConcentration")
UMI_Glo2 = subset(UMI_Glo2, subset =  nCount_RNA > 200)

UMI_Glo2 = NormalizeData(UMI_Glo2, normalization.method = "RC")
UMI_Glo2@assays$RNA@data = as.matrix(log2(UMI_Glo2@assays$RNA@data + 1))
UMI_Glo = CreateSeuratObject(counts = UMI_Glo2@assays$RNA@data)
UMI_Glo@assays$RNA@counts = UMI_Glo2@assays$RNA@counts
UMI_Glo@assays$RNA@data = UMI_Glo2@assays$RNA@data
UMI_Glo@meta.data$mRNAConcentration = UMI_Glo2@meta.data$mRNAConcentration

Idents(object = UMI_Glo) = UMI_Glo@meta.data$mRNAConcentration

#Panel A/C Top, 600 x 350
FeatureScatter(
  UMI_Glo,
  'nCount_RNA',
  'mRNAConcentration',
  slot = "data",
  pt.size = 1.5,
  cols = scattercol
) + xlab("Transcripts") + ylab("mRNA (pg)") +
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

UMI_Glo_curve = as.data.frame(t(as.data.frame(UMI_Glo@assays$RNA@data)))
UMI_Glo_curve[, "celltype"] = UMI_Glo@active.ident
UMI_Glo_curve = aggregate(. ~ celltype, UMI_Glo_curve, mean)
rownames(UMI_Glo_curve) = UMI_Glo_curve$celltype
UMI_Glo_curve$celltype = NULL
UMI_Glo_curve = as.data.frame(t(UMI_Glo_curve))
UMI_Glo_curve[, "total"] = rowSums(UMI_Glo_curve)
UMI_Glo_curve[, "genes"] = rownames(UMI_Glo_curve)
UMI_Glo_curve = UMI_Glo_curve[order(UMI_Glo_curve$total), ]

order.total = order(UMI_Glo_curve$total, decreasing = TRUE)
UMI_Glo_curve$rank <- NA
UMI_Glo_curve$rank[order.total] <- 1:nrow(UMI_Glo_curve)

#Panel B/D Top, 600 x 350
UMI_Glo_curve2 = UMI_Glo_curve[UMI_Glo_curve$rank <= 20000, ]
UMI_Glo_tidy = UMI_Glo_curve2 %>%
  select(`3.75`, `7.5`, `15`, `30`, rank) %>%
  gather(key = "variable", value = "value",-rank)

ggplot(UMI_Glo_tidy, aes(x = rank, y = value)) +
  geom_smooth(
    aes(color = variable),
    method = "gam",
    size = 1.5,
    alpha = 1
  ) +
  theme_classic() +
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
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = c("seagreen2", "seagreen4", "seagreen1", "seagreen3")) +
  NoLegend()


#UMI Absolute scaling by raw UMIs
UMI_Abs <- CreateSeuratObject(counts = umi_raw)
UMI_Abs = AddMetaData(UMI_Abs,
                      metadata = barcodes$mRNA_amount,
                      col.name = "mRNAConcentration")
UMI_Abs = subset(UMI_Abs, subset =  nCount_RNA > 200)

UMI_Abs@assays$RNA@data = as.matrix(log2(UMI_Abs@assays$RNA@counts + 1))
Idents(object = UMI_Abs) = UMI_Abs@meta.data$mRNAConcentration

#Panel A/C Middle, 600 x 350
FeatureScatter(
  UMI_Abs,
  'nCount_RNA',
  'mRNAConcentration',
  slot = "data",
  pt.size = 1.5,
  cols = scattercol
) + xlab("Transcripts") + ylab("mRNA (pg)") +
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

UMI_Abs_curve = as.data.frame(t(as.data.frame(UMI_Abs@assays$RNA@data)))
UMI_Abs_curve[, "celltype"] = UMI_Abs@active.ident
UMI_Abs_curve = aggregate(. ~ celltype, UMI_Abs_curve, mean)
rownames(UMI_Abs_curve) = UMI_Abs_curve$celltype
UMI_Abs_curve$celltype = NULL
UMI_Abs_curve = as.data.frame(t(UMI_Abs_curve))
UMI_Abs_curve[, "total"] = rowSums(UMI_Abs_curve)
UMI_Abs_curve[, "genes"] = rownames(UMI_Abs_curve)
UMI_Abs_curve = UMI_Abs_curve[order(UMI_Abs_curve$total), ]

order.total = order(UMI_Abs_curve$total, decreasing = TRUE)
UMI_Abs_curve$rank <- NA
UMI_Abs_curve$rank[order.total] <- 1:nrow(UMI_Abs_curve)

UMI_Abs_curve2 = UMI_Abs_curve
UMI_Glo_curve2 = UMI_Glo_curve2[UMI_Glo_curve2$rank <= 20000, ]
UMI_Abs_curve2$`3.75` = log2(UMI_Abs_curve2$`3.75` + 1)
UMI_Abs_curve2$`7.5` = log2(UMI_Abs_curve2$`7.5` + 1)
UMI_Abs_curve2$`15` = log2(UMI_Abs_curve2$`15` + 1)
UMI_Abs_curve2$`30` = log2(UMI_Abs_curve2$`30` + 1)

UMI_Abs_tidy = UMI_Abs_curve2 %>%
  select(`3.75`, `7.5`, `15`, `30`, rank) %>%
  gather(key = "variable", value = "value",-rank)

#Panel B/D Middle, 600 x 350
ggplot(UMI_Abs_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable),
              size = 1.5,
              alpha = 1,
              se = FALSE) +
  theme_classic() +
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
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = c("seagreen2", "seagreen4", "seagreen1", "seagreen3")) +
  NoLegend()

#ERCC Absolute scaling by scran/ERCC
ERCC_Abs2 <- CreateSeuratObject(counts = ercc_raw)
ERCC_Abs2 = AddMetaData(ERCC_Abs2,
                        metadata = barcodes$mRNA_amount,
                        col.name = "mRNAConcentration")
ERCC_Abs2 = AddMetaData(ERCC_Abs2,
                        PercentageFeatureSet(ERCC_Abs2, pattern = "^ERCC-"),
                        col.name = "ERCC")
ERCC_Abs2 = subset(ERCC_Abs2, subset = nCount_RNA > 200 &
                     ERCC > 0.1)

filteredcounts = GetAssayData(object = ERCC_Abs2,
                              assay = "RNA",
                              slot = "data")
scran.data = SingleCellExperiment(assays = list(counts = filteredcounts))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <-
  splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))

scran.data <- computeSpikeFactors(scran.data, "ERCC")
scran.data = logNormCounts(scran.data, log = FALSE)

ERCC_Abs <-
  CreateSeuratObject(counts = as.sparse(x = assay(scran.data, "normcounts")))
ERCC_Abs@meta.data$mRNAConcentration = ERCC_Abs2@meta.data$mRNAConcentration
ERCC_Abs = AddMetaData(ERCC_Abs,
                       PercentageFeatureSet(ERCC_Abs, pattern = "^ERCC-"),
                       col.name = "ERCC")

ERCC_Abs@assays$RNA@data = as.matrix(log2(ERCC_Abs@assays$RNA@counts + 1))
Idents(object = ERCC_Abs) = ERCC_Abs@meta.data$mRNAConcentration

#Panel A/C Bottom, 600 x 350
FeatureScatter(
  ERCC_Abs,
  'nCount_RNA',
  'mRNAConcentration',
  slot = "data",
  pt.size = 1.5,
  cols = scattercol
) + xlab("Transcripts") + ylab("mRNA (pg)") +
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

ERCC_Abs_curve = as.data.frame(t(as.data.frame(ERCC_Abs@assays$RNA@data)))
ERCC_Abs_curve[, "celltype"] = ERCC_Abs@active.ident
ERCC_Abs_curve = aggregate(. ~ celltype, ERCC_Abs_curve, mean)
rownames(ERCC_Abs_curve) = ERCC_Abs_curve$celltype
ERCC_Abs_curve$celltype = NULL
ERCC_Abs_curve = as.data.frame(t(ERCC_Abs_curve))
ERCC_Abs_curve[, "total"] = rowSums(ERCC_Abs_curve)
ERCC_Abs_curve[, "genes"] = rownames(ERCC_Abs_curve)
ERCC_Abs_curve = ERCC_Abs_curve[order(ERCC_Abs_curve$total), ]

order.total = order(ERCC_Abs_curve$total, decreasing = TRUE)
ERCC_Abs_curve$rank <- NA
ERCC_Abs_curve$rank[order.total] <- 1:nrow(ERCC_Abs_curve)

ERCC_Abs_curve2 = ERCC_Abs_curve[ERCC_Abs_curve$rank <= 20000, ]
ERCC_Abs_tidy = ERCC_Abs_curve2 %>%
  select(`3.75`, `7.5`, `15`, `30`, rank) %>%
  gather(key = "variable", value = "value",-rank)

#Panel B/D Bottom, 600 x 350
ggplot(ERCC_Abs_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable),
              size = 1.5,
              alpha = 1,
              se = FALSE) +
  theme_classic() +
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
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = c("seagreen2", "seagreen4", "seagreen1", "seagreen3")) +
  NoLegend()

#UMI/ERCC Comparison
UMI_Abs_transcripts = data.frame(
  "UMI Transcripts" = UMI_Abs@meta.data$nCount_RNA,
  "Names" = colnames(UMI_Abs@assays$RNA@data),
  "Value" = UMI_Abs@meta.data$mRNAConcentration
)
ERCC_Abs_transcripts = data.frame(
  "ERCC Transcripts" = ERCC_Abs@meta.data$nCount_RNA,
  "Names" = colnames(ERCC_Abs@assays$RNA@data),
  "Value" = ERCC_Abs@meta.data$mRNAConcentration
)
comp_Abs = merge(UMI_Abs_transcripts, ERCC_Abs_transcripts, by = "Names")
comp_Abs$UMI.Transcripts = normalize(comp_Abs$UMI.Transcripts,
                                     method = "range",
                                     range = c(0, 10))
comp_Abs$ERCC.Transcripts = normalize(comp_Abs$ERCC.Transcripts,
                                      method = "range",
                                      range = c(0, 10))

#Panel E/F Bottom, 600 x 400
library(ggpubr)
ggplot(comp_Abs, aes(x = ERCC.Transcripts, y = UMI.Transcripts)) +
  geom_point(size = 1.5, aes(color = as.factor(Value.x))) +
  geom_smooth(
    color = "black",
    fill = "grey90",
    linetype = "dashed",
    method = "lm"
  ) +
  theme_classic() +
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
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "UMI-Normalization", x = "ERCC-Normalization") +
  scale_color_manual(values = c("seagreen4", "seagreen3", "seagreen2", "seagreen1")) +
  xlim(0, 10) +
  ylim(0, 10) +
  stat_regline_equation(size = 4.5) +
  NoLegend() #600 x 400
