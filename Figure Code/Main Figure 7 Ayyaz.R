library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(scran)
library(scRNAseq)
library(scater)
library(DoubletFinder)

crypt <- readRDS("Ayyaz2019/crypt.rds") # UMI-normalized with combined WT and IR from original count matrices
crypt@assays$RNA@data <- as.matrix(log2(crypt@assays$RNA@counts + 1))
crypt <- FindVariableFeatures(crypt)
crypt <- ScaleData(crypt)
crypt <- RunPCA(crypt)
crypt <- FindNeighbors(crypt, dims = 1:40)
crypt <- FindClusters(crypt, resolution = 0.4)
crypt <- RunUMAP(crypt, dims = 1:40, min.dist = 0.75)
crypt@meta.data$logRNA <- log2(crypt@meta.data$nCount_RNA + 1)

# Figure Panels
new.cluster.ids <-
  c(
    "Enterocyte-1",
    "Lgr5+ Stem Cell",
    "Enterocyte-2",
    "Enterocyte-3",
    "Lymphocyte-1",
    "Enterocyte-4",
    "Lymphocyte-2",
    "Enterocyte-5",
    "Goblet Cell-1",
    "Goblet Cell-2",
    "Clu+ Stem Cell",
    "Crypt-base Columnar-2",
    "Enterocyte-6",
    "Entero-Endocrine",
    "Tuft Cell",
    "Paneth Cell"
  )

names(new.cluster.ids) <- levels(crypt)
crypt <- RenameIdents(crypt, new.cluster.ids)
irradiated <- subset(x = crypt, Experiment == "Irradiated")
normal <- subset(x = crypt, Experiment == "Normal")

# Panel F, 600 x 350
FeaturePlot(normal, features = c("logRNA"), pt.size = 0.5) +
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
  ) + labs(title = "log2 UMI")
FeaturePlot(normal, features = c("Clu"), pt.size = 0.5, ) +
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
  ) + labs(title = "Clu")
FeaturePlot(normal, features = c("Lgr5"), pt.size = 0.5) +
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
  ) + labs(title = "Lgr5")

# Panel G, 600 x 350
FeaturePlot(irradiated,
  features = c("logRNA"),
  pt.size = 0.5
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
    )
  ) + labs(title = "log2 UMI")
FeaturePlot(irradiated, features = c("Clu"), pt.size = 0.5) +
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
  ) + labs(title = "Clu")
FeaturePlot(irradiated, features = c("Lgr5"), pt.size = 0.5) +
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
  ) + labs(title = "Lgr5")

ordering <- data.frame(RNA = crypt$logRNA, celltype = crypt@active.ident)
ordering <- aggregate(. ~ celltype, ordering, median)
ordering <- ordering[order(-ordering$RNA), ]

Idents(crypt) <- factor(Idents(crypt), levels = ordering$celltype)

# Panel H Top, 600 x 350
organcols <- c("grey", "orangered2", rep("grey", 14))
VlnPlot(
  crypt,
  features = c("logRNA"),
  cols = organcols,
  ncol = 1,
  pt.size = 0
) +
  ylab("log2 UMIs") + xlab("") + NoLegend() +
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 10,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 0.5,
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
    axis.text.x = element_text(angle = 45)
  )

# Cell Cycle Scoring
library(stringr)
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
crypt <- NormalizeData(crypt, normalization.method = "RC")
crypt@assays$RNA@data <- as.matrix(log2(crypt@assays$RNA@data + 1))
crypt <-
  CellCycleScoring(
    crypt,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE
  )

# Panel H Middle, 600 x 350
dimcols <- rep("seagreen3", 16)
VlnPlot(
  crypt,
  features = c("S.Score"),
  cols = dimcols,
  ncol = 1,
  pt.size = 0
) +
  ylab("S-Phase Score") + xlab("") + NoLegend() +
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 10,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 0.5,
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
    axis.text.x = element_text(angle = 45)
  )

# Panel H Bottom, 600 x 350
dimcols <- rep("lightblue", 16)
VlnPlot(
  crypt,
  features = c("G2M.Score"),
  cols = dimcols,
  ncol = 1,
  pt.size = 0
) +
  ylab("G2/M-Phase Score") + xlab("") + NoLegend() +
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 10,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 0.5,
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
    axis.text.x = element_text(angle = 45)
  )
