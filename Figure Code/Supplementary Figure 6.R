library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)

tabulafacs <- readRDS("Tabula FACS/tabulafacs.rds")
tabula_unnorm <- readRDS("Tabula FACS/tabulaunnorm.rds")
tabula_unnorm@meta.data$facsRNA <- tabulafacs@meta.data$nCount_RNA
cols <- read.csv(file = "Tabula FACS/tabulafacscols.csv") # consistent colors across figures

autoprocess <- function(subset_organ) {
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

# Find shared cell types
unique_cellpairs <- as.data.frame(tabula_unnorm@meta.data$celltype)
colnames(unique_cellpairs) <- "celltype"
unique_cellpairs$tissue <- as.data.frame(tabula_unnorm@meta.data$tissue)
unique_cellpairs <- unique_cellpairs[!duplicated(t(apply(unique_cellpairs, 1, sort))), ]

# T-cells
t_cell <- subset(tabula_unnorm, subset = celltype == "T cell")
t_cell <- autoprocess(t_cell)

DimPlot(
  t_cell,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Fat" = "#7CAE00",
    "Limb_Muscle" = "#00C1A3",
    "Lung" = "#00BAE0",
    "Spleen" = "#E76BF3"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "T cells")

t_cell_df <- data.frame(
  tissue = t_cell@meta.data$tissue,
  transcripts = t_cell@meta.data$facsRNA
)
t_cell_df$colors <- cols$colors[match(t_cell_df$tissue, cols$Organs)]
t_cell_df$tissue <- gsub("_", " ", t_cell_df$tissue)

ggplot(t_cell_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts
)) +
  theme_classic() +
  geom_boxplot(aes(fill = colors)) +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  scale_fill_identity() +
  labs(title = "", y = "Transcripts", x = "") +
  NoLegend()


# B-Cells
b_cell <- subset(tabula_unnorm, subset = celltype == "B cell")
b_cell <- autoprocess(b_cell)

DimPlot(
  b_cell,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Fat" = "#7CAE00",
    "Limb_Muscle" = "#00C1A3",
    "Lung" = "#00BAE0",
    "Spleen" = "#E76BF3",
    "Marrow" = "#35A2FF",
    "Liver" = "#00BFC4"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "B cells")

b_cell_df <- data.frame(
  tissue = b_cell@meta.data$tissue,
  transcripts = b_cell@meta.data$facsRNA
)
b_cell_df$colors <- cols$colors[match(b_cell_df$tissue, cols$Organs)]
b_cell_df$tissue <- gsub("_", " ", b_cell_df$tissue)

ggplot(b_cell_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts
)) +
  theme_classic() +
  geom_boxplot(aes(fill = colors)) +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  scale_fill_identity() +
  labs(title = "", y = "Transcripts", x = "") +
  NoLegend()

# Macrophages
macrophage <- subset(tabula_unnorm, subset = celltype == "macrophage")
macrophage <- autoprocess(macrophage)

DimPlot(
  macrophage,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Brain_Myeloid" = "#D89000",
    "Limb_Muscle" = "#00C1A3",
    "Diaphragm" = "#A3A500",
    "Spleen" = "#E76BF3",
    "Marrow" = "#35A2FF",
    "Kidney" = "#00BB4E"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "Macrophages")

macrophage_df <- data.frame(
  tissue = macrophage@meta.data$tissue,
  transcripts = macrophage@meta.data$facsRNA
)
macrophage_df$organs <- cols$colors[match(macrophage_df$tissue, cols$Organs)]
macrophage_df$tissue <- gsub("_", " ", macrophage_df$tissue)

ggplot(macrophage_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts,
  fill = organs
)) +
  geom_boxplot() +
  theme_classic() +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  labs(title = "", y = "Transcripts", x = "") +
  scale_fill_identity() +
  NoLegend()

# Endothelial Cells
endo <- subset(tabula_unnorm, subset = celltype == "endothelial cell")
endo <- autoprocess(endo)

DimPlot(
  endo,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Pancreas" = "#9590FF",
    "Brain_Non-Myeloid" = "#D89000",
    "Mammary Gland" = "#00B0F6",
    "Trachea" = "#FF6A98",
    "Aorta" = "#F8766D",
    "Fat" = "#7CAE00",
    "Kidney" = "#00BB4E",
    "Heart" = "#39B600",
    "Limb_Muscle" = "#00C1A3",
    "Diaphragm" = "#A3A500"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "Endothelial Cells")

endo_df <- data.frame(
  tissue = endo@meta.data$tissue,
  transcripts = endo@meta.data$facsRNA
)
endo_df$colors <- cols$colors[match(endo_df$tissue, cols$Organs)]
endo_df$tissue <- gsub("_", " ", endo_df$tissue)

ggplot(endo_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts
)) +
  theme_classic() +
  geom_boxplot(aes(fill = colors)) +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  labs(title = "", y = "Transcripts", x = "") +
  scale_fill_identity() +
  NoLegend()

# Droplet Datasets
tabula10X <- readRDS("Tabula 10X/tabula10X.rds")
unique_cellpairs <- as.data.frame(tabula10X@meta.data$celltype)
colnames(unique_cellpairs) <- "celltype"
unique_cellpairs$tissue <- as.data.frame(tabula10X@meta.data$tissue)
unique_cellpairs <- unique_cellpairs[!duplicated(t(apply(unique_cellpairs, 1, sort))), ]
cols10X <- read.csv("Tabula 10X/tabula10Xcols.csv")

# T-cells
t_cell <- subset(tabula10X, subset = celltype == "T cell")
t_cell <- autoprocess(t_cell)

DimPlot(
  t_cell,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Marrow" = "#00B4F0",
    "Mammary Gland" = "#00BFC4",
    "Lung" = "#00C08B",
    "Limb Muscle" = "#7CAE00",
    "Spleen" = "#619CFF"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "T cells")

t_cell_df <- data.frame(
  tissue = t_cell@meta.data$tissue,
  transcripts = t_cell@meta.data$nCount_RNA
)
t_cell_df$colors <- cols10X$colors[match(t_cell_df$tissue, cols10X$Organs)]
t_cell_df$tissue <- gsub("_", " ", t_cell_df$tissue)

ggplot(t_cell_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts
)) +
  theme_classic() +
  geom_boxplot(aes(fill = colors)) +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  scale_fill_identity() +
  labs(title = "", y = "Transcripts", x = "") +
  NoLegend()

# B-Cells
b_cell <- subset(tabula10X, subset = celltype == "B cell")
b_cell <- autoprocess(b_cell)

DimPlot(
  b_cell,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Mammary Gland" = "#00BFC4",
    "Spleen" = "#619CFF",
    "Lung" = "#00C08B",
    "Limb Muscle" = "#7CAE00"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "B cells")

b_cell_df <- data.frame(
  tissue = b_cell@meta.data$tissue,
  transcripts = b_cell@meta.data$nCount_RNA
)
b_cell_df$colors <- cols10X$colors[match(b_cell_df$tissue, cols10X$Organs)]
b_cell_df$tissue <- gsub("_", " ", b_cell_df$tissue)

ggplot(b_cell_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts
)) +
  theme_classic() +
  geom_boxplot(aes(fill = colors)) +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  scale_fill_identity() +
  labs(title = "", y = "Transcripts", x = "") +
  NoLegend()

# Macrophage
macrophage <- subset(tabula10X, subset = celltype == "macrophage")
macrophage <- autoprocess(macrophage)
DimPlot(
  macrophage,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Mammary Gland" = "#00BFC4",
    "Spleen" = "#619CFF",
    "Limb Muscle" = "#7CAE00",
    "Kidney" = "#B79F00",
    "Marrow" = "#00B4F0"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "Macrophages")

macrophage_df <- data.frame(
  tissue = macrophage@meta.data$tissue,
  transcripts = macrophage@meta.data$nCount_RNA
)
macrophage_df$organs <- cols10X$colors[match(macrophage_df$tissue, cols10X$Organs)]
macrophage_df$tissue <- gsub("_", " ", macrophage_df$tissue)

ggplot(macrophage_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts,
  fill = organs
)) +
  geom_boxplot() +
  theme_classic() +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  labs(title = "", y = "Transcripts", x = "") +
  scale_fill_identity() +
  NoLegend()

# Endothelial Cell
endo <- subset(tabula10X, subset = celltype == "endothelial cell")
endo <- autoprocess(endo)
DimPlot(
  endo,
  reduction = "umap",
  pt.size = 1,
  group.by = "tissue",
  cols = c(
    "Mammary Gland" = "#00BFC4",
    "Bladder" = "#F8766D",
    "Heart and Aorta" = "#DE8C00",
    "Trachea" = "#FF64B0",
    "Limb Muscle" = "#7CAE00"
  )
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) +
  labs(title = "Endothelial Cells")

endo_df <- data.frame(
  tissue = endo@meta.data$tissue,
  transcripts = endo@meta.data$nCount_RNA
)
endo_df$colors <- cols10X$colors[match(endo_df$tissue, cols10X$Organs)]
endo_df$tissue <- gsub("_", " ", endo_df$tissue)

ggplot(endo_df, aes(
  x = reorder(tissue, -transcripts, .fun = "median"),
  y = transcripts
)) +
  theme_classic() +
  geom_boxplot(aes(fill = colors)) +
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
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5)
  ) +
  labs(title = "", y = "Transcripts", x = "") +
  scale_fill_identity() +
  NoLegend()
