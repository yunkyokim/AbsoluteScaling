library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)

tabulafacs <- readRDS("Tabula FACS/tabulafacs.rds")
tabula_unnorm <- readRDS("Tabula FACS/tabulaunnorm.rds")
tabula_unnorm@meta.data$facsRNA <- tabulafacs@meta.data$nCount_RNA

autoprocessGLO <- function(subset_organ) {
  subset_organ <- NormalizeData(subset_organ, normalization.method = "RC")
  subset_organ@assays$RNA@data <- as.matrix(log2(subset_organ@assays$RNA@data +
    1))
  subset_organ <-
    FindVariableFeatures(subset_organ,
      selection.method = "vst",
      nfeatures = 2000
    )
  all.genes <- rownames(subset_organ)
  subset_organ <- ScaleData(subset_organ, features = all.genes)
  subset_organ <-
    RunPCA(subset_organ, features = VariableFeatures(object = subset_organ))
  subset_organ <- FindNeighbors(subset_organ, dims = 1:15)
  subset_organ <- FindClusters(subset_organ, resolution = 0.4)
  subset_organ <-
    RunUMAP(subset_organ, dims = 1:15, min.dist = 0.75)
  return(subset_organ)
}
autoprocessABS <- function(subset_organ) {
  subset_organ@assays$RNA@data <- as.matrix(log2(subset_organ@assays$RNA@counts +
    1))
  subset_organ <-
    FindVariableFeatures(subset_organ,
      selection.method = "vst",
      nfeatures = 2000
    )
  all.genes <- rownames(subset_organ)
  subset_organ <- ScaleData(subset_organ, features = all.genes)
  subset_organ <-
    RunPCA(subset_organ, features = VariableFeatures(object = subset_organ))
  subset_organ <- FindNeighbors(subset_organ, dims = 1:15)
  subset_organ <- FindClusters(subset_organ, resolution = 0.4)
  subset_organ <-
    RunUMAP(subset_organ, dims = 1:15, min.dist = 0.75)
  return(subset_organ)
}

# Hepatocyte Polyploidy
tabula_liver <- subset(tabula_unnorm, subset = celltype == "hepatocyte")
tabula_liver <- autoprocessGLO(tabula_liver)

# Panel A
FeaturePlot(tabula_liver,
  features = c("facsRNA"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Transcripts (Absolute Scaling)")

# Panel B
FeaturePlot(tabula_liver,
  features = c("Mlxipl"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Mlxipl (Global Scaling)") # 650 x 350

# Panel C
FeaturePlot(tabula_liver,
  features = c("Nr1i3"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Nr1i3 (Global Scaling)") # 650 x 350

# Panel D
FeaturePlot(tabula_liver,
  features = c("Lifr"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Lifr (Global Scaling)") # 650 x 350

# Bladder Polyploidy
tabula_bladder <- subset(tabula_unnorm, subset = celltype == "bladder urothelial cell")
tabula_bladder <- autoprocessGLO(tabula_bladder)

# Panel F
DimPlot(tabula_bladder, reduction = "umap", pt.size = 1) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "Bladder Epithelium")

# Panel E
FeaturePlot(tabula_bladder,
  features = c("facsRNA"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Transcripts (Absolute Scaling)")

# Panel H
FeaturePlot(tabula_bladder,
  features = c("Krt20"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Krt20 (Global Scaling)")

# Panel J
FeaturePlot(tabula_bladder,
  features = c("Trp63"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Trp63 (Global Scaling)")

# Panel L
FeaturePlot(tabula_bladder,
  features = c("Upk2"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Upk2 (Global Scaling)") # 650 x 350

# Panel N
FeaturePlot(tabula_bladder,
  features = c("Krt14"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Krt14 (Global Scaling)")

tabula_bladder <- RenameIdents(
  object = tabula_bladder,
  `0` = "Superficial",
  `1` = "Intermediate 1",
  `2` = "Intermediate 2",
  `3` = "Basal"
)
tabula_bladder@meta.data$boxplotclusters <- tabula_bladder@active.ident
auto_boxplot <- function(organ, gene, title, col) {
  organboxplot <- data.frame(
    features = as.data.frame(FetchData(organ, gene, slot = "data")),
    clusters = organ@meta.data$boxplotclusters
  )
  colnames(organboxplot) <- c("features", "clusters")
  ggplot(organboxplot, aes(x = clusters, y = features)) +
    theme_classic() +
    geom_boxplot(fill = col) +
    theme(
      aspect.ratio = 1,
      text = element_text(size = 12),
      plot.margin = margin(
        t = 7,
        r = 0,
        b = 7,
        l = 0,
        unit = "pt"
      ),
      axis.title = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 10, colour = "black"),
      plot.title = element_text(size = 11, hjust = 0.5)
    ) +
    labs(title = title, y = "Expression", x = "") +
    NoLegend()
}

# Panels I/K/M/O
auto_boxplot(tabula_bladder, "facsRNA", "Transcripts", "thistle3")
auto_boxplot(tabula_bladder, "Krt20", "Krt20", "thistle3")
auto_boxplot(tabula_bladder, "Upk2", "Upk2", "thistle3")
auto_boxplot(tabula_bladder, "Trp63", "Trp63", "thistle3")
auto_boxplot(tabula_bladder, "Krt14", "Krt14", "thistle3")

# Panel G
FeaturePlot(tabula_bladder,
  features = c("nCount_RNA"),
  pt.size = 1
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Transcripts")
