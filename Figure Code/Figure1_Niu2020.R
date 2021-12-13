library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(viridis)
library(cowplot)

#setwd("current_dir/")

#Data import, files available upon request
#Import Seurat object
gonads = readRDS("gonadstotal.rds")
gonads@meta.data$noclusters = 0
gonads = subset(gonads, subset = Car2 == 0 &
                  Cldn5 == 0) #Immune cell removal

#Standard Normalization/Processing
gonads@assays$RNA@data = as.matrix(log2(gonads@assays$RNA@counts + 1))
gonads <- FindVariableFeatures(gonads)
all.genes <- rownames(gonads)
gonads <- ScaleData(gonads, features = all.genes)
gonads <-
  RunPCA(
    gonads,
    npcs = 100,
    ndims.print = 1:5,
    nfeatures.print = 5
  )
gonads <- FindNeighbors(gonads, dims = 1:15)
gonads <- FindClusters(gonads, resolution = 2)
gonads <- RunUMAP(gonads, dims = 1:15, min.dist = 0.75)

#PGC Assignment
Oct4 = as.data.frame(gonads@assays$RNA@counts["Pou5f1", ])
colnames(Oct4) = "Oct4"
for (row in 1:nrow(Oct4)) {
  expr <- Oct4[row, "Oct4"]
  if (expr > 0.5) {
    Oct4[row, "status"] = "Oct4+"
  } else {
    Oct4[row, "status"] = "Oct4-"
  }
}
gonads@meta.data$germ = Oct4$status

#Figure Panels
dimcols = c("grey80", "violetred2")
dimcols2 = rep(dimcols, 4)

gonads_log = NormalizeData(gonads, normalization.method = "RC")
gonads_log = CreateSeuratObject(counts = gonads_log@assays$RNA@data)
gonads_log@assays$RNA@data = as.matrix(log2(gonads_log@assays$RNA@counts + 1))
gonads_log@meta.data$nCount_logRNA = colSums(gonads_log@assays$RNA@data)
gonads_log@meta.data$germ = gonads@meta.data$germ
gonads_log@meta.data$timepoint = gonads@meta.data$timepoint

#Panel O, 700 x 300
VlnPlot(
  gonads_log,
  features = c("nCount_logRNA"),
  split.plot = TRUE,
  group.by = "timepoint",
  cols = dimcols,
  split.by = "germ",
  pt.size = 0
) +
  scale_y_continuous(
    labels = function(x)
      format(x, scientific = TRUE)
  ) +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 10,
    colour = dimcols2,
    shape = 95
  ) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35,
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
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
#Panel P, 700 x 300
VlnPlot(
  gonads,
  features = c("nCount_RNA"),
  split.plot = TRUE,
  group.by = "timepoint",
  cols = dimcols,
  split.by = "germ",
  pt.size = 0
) +
  scale_y_continuous(
    labels = function(x)
      format(x, scientific = TRUE)
  ) +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 10,
    colour = dimcols2,
    shape = 95
  ) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35,
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
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#Cell Cycle Scoring
library(stringr)
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
gonads_log <-
  CellCycleScoring(
    gonads_log,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE
  )

#Panel Q Top, 700 x 300
dimcols = c("grey80", "seagreen3")
dimcols2 = rep(dimcols, 4)
VlnPlot(
  gonads_log,
  features = c("S.Score"),
  split.plot = TRUE,
  group.by = "timepoint",
  cols = dimcols,
  split.by = "germ",
  pt.size = 0
) +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 10,
    colour = dimcols2,
    shape = 95
  ) +
  ylab("S-Phase Score") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35,
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
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#Panel Q Bottom, 700 x 300
dimcols = c("grey80", "lightblue")
dimcols2 = rep(dimcols, 4)
VlnPlot(
  gonads_log,
  features = c("G2M.Score"),
  split.plot = TRUE,
  group.by = "timepoint",
  cols = dimcols,
  split.by = "germ",
  pt.size = 0
) +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 10,
    colour = dimcols2,
    shape = 95
  ) +
  ylab("G2/M-Phase Score") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 0.35,
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
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


#Subset of E12.5
gonad_E12 = subset(gonads, subset = timepoint == "E12.5")
gonad_E12 <- FindVariableFeatures(gonad_E12)
all.genes <- rownames(gonad_E12)
gonad_E12 <- ScaleData(gonad_E12, features = all.genes)
gonad_E12 <-
  RunPCA(
    gonad_E12,
    npcs = 100,
    ndims.print = 1:5,
    nfeatures.print = 5
  )
gonad_E12 <- FindNeighbors(gonad_E12, dims = 1:15)
gonad_E12 <- FindClusters(gonad_E12, resolution = 0.5)
gonad_E12 <- RunUMAP(gonad_E12, dims = 1:15, min.dist = 0.75)

Oct4 = as.data.frame(gonad_E12@assays$RNA@counts["Pou5f1", ])
colnames(Oct4) = "Oct4"
for (row in 1:nrow(Oct4)) {
  expr <- Oct4[row, "Oct4"]
  
  if (expr > 0.3) {
    Oct4[row, "status"] = "Oct4+"
  } else {
    Oct4[row, "status"] = "Oct4-"
  }
}
gonad_E12@meta.data$germ = Oct4$status

#Fold-change calculation
gonad_E12_transcripts = data.frame(transcripts = gonad_E12@meta.data$nCount_RNA,
                                   identity = gonad_E12@meta.data$germ)
gonad_E12_transcripts = aggregate(. ~ identity, gonad_E12_transcripts, median)
rownames(gonad_E12_transcripts) = gonad_E12_transcripts$identity
gonad_E12_transcripts["Oct4+", "transcripts"] / gonad_E12_transcripts["Oct4-", "transcripts"] #fold-change value, = 1.622191

gonad_E12_log = NormalizeData(gonad_E12, normalization.method = "RC")
gonad_E12_log = CreateSeuratObject(counts = gonad_E12_log@assays$RNA@data)
gonad_E12_log@assays$RNA@data = as.matrix(log2(gonad_E12_log@assays$RNA@counts + 1))
gonad_E12_log@meta.data$germ = gonad_E12@meta.data$germ
gonad_E12_log <- FindVariableFeatures(gonad_E12_log)
all.genes <- rownames(gonad_E12_log)
gonad_E12_log <- ScaleData(gonad_E12_log, features = all.genes)
gonad_E12_log <-
  RunPCA(
    gonad_E12_log,
    npcs = 100,
    ndims.print = 1:5,
    nfeatures.print = 5
  )

gonad_E12_log <- FindNeighbors(gonad_E12_log, dims = 1:15)
gonad_E12_log <- FindClusters(gonad_E12_log, resolution = 0.5)
gonad_E12_log <-
  RunUMAP(gonad_E12_log, dims = 1:15, min.dist = 0.75)
gonad_E12_log@meta.data$nCount_logRNA = colSums(gonad_E12_log@assays$RNA@data)

#Figure Panels
#Panel K, 600 x 350
DimPlot(
  gonad_E12_log,
  reduction = "umap",
  group.by = 'germ',
  cols = dimcols,
  pt.size = 0.1
)  +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    plot.title = element_blank(),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) + NoLegend()

#Panel L, 800 x 300
DimPlot(
  gonad_E12,
  reduction = "umap",
  group.by = 'germ',
  cols = dimcols,
  pt.size = 0.1
)  +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    plot.title = element_blank(),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) + NoLegend()

#Panel M, 600 x 350
FeaturePlot(gonad_E12_log,
            features = c("nCount_logRNA"),
            pt.size = 0.1)  +
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

#Panel N, 600 x 350
FeaturePlot(gonad_E12,
            features = c("nCount_RNA"),
            pt.size = 0.1)  +
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