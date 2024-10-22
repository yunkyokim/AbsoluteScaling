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
library(monocle3)

marrow <- readRDS("Tabula 10X BM/marrow.rds") # subset from Tabula 10X/tabula10X.rds
marrow@assays$RNA@data <- as.matrix(log2(marrow@assays$RNA@counts + 1))
marrow <-
  FindVariableFeatures(marrow, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(marrow)
marrow <- ScaleData(marrow, features = all.genes)
marrow <-
  RunPCA(marrow, features = VariableFeatures(object = marrow), npcs = 150)
marrow <- RunUMAP(marrow, dims = 1:100, min.dist = 0.75)
marrow <- FindNeighbors(marrow, dims = 1:100)
marrow <- FindClusters(marrow)

# Panel A
DimPlot(marrow,
  reduction = "umap",
  pt.size = 0.1,
  group.by = "celltype"
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Marrow Tabula 10X")

# Panel B
FeaturePlot(marrow,
  features = c("nCount_RNA"),
  pt.size = 0.05
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Transcripts")

# Subset monocytes
monocyte <-
  subset(
    marrow,
    subset = celltype == "hematopoietic precursor cell" |
      celltype == "promonocyte" |
      celltype == "monocyte"
  )
monocyte@assays$RNA@data <- as.matrix(log2(monocyte@assays$RNA@counts + 1))
monocyte <- NormalizeData(monocyte)
monocyte <-
  FindVariableFeatures(monocyte, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(monocyte)
monocyte <- ScaleData(monocyte, features = all.genes)
monocyte <-
  RunPCA(monocyte,
    features = VariableFeatures(object = monocyte),
    npcs = 150
  )
monocyte <- RunUMAP(monocyte, dims = 1:60, min.dist = 0.4)
monocyte <- FindNeighbors(monocyte, dims = 1:60)
monocyte <- FindClusters(monocyte)

monocyte <- subset(monocyte, subset = seurat_clusters != 8)
monocyte <- subset(monocyte, subset = seurat_clusters != 10)
monocyte <- subset(monocyte, subset = seurat_clusters != 2)
monocyte <- subset(monocyte, subset = seurat_clusters != 3)
monocyte <-
  RunPCA(monocyte,
    features = VariableFeatures(object = monocyte),
    npcs = 150
  )
monocyte <- RunUMAP(monocyte, dims = 1:60, min.dist = 0.4)
monocyte <- FindNeighbors(monocyte, dims = 1:60)
monocyte <- FindClusters(monocyte)

# Subset Erythroblasts
erythroblast <-
  subset(
    marrow,
    subset = celltype == "hematopoietic precursor cell" |
      celltype == "proerythroblast" |
      celltype == "erythroblast"
  )
erythroblast@assays$RNA@data <- as.matrix(log2(erythroblast@assays$RNA@counts +
  1))
erythroblast <- NormalizeData(erythroblast)
erythroblast <-
  FindVariableFeatures(erythroblast,
    selection.method = "vst",
    nfeatures = 5000
  )
all.genes <- rownames(erythroblast)
erythroblast <- ScaleData(erythroblast, features = all.genes)
erythroblast <-
  RunPCA(erythroblast,
    features = VariableFeatures(object = erythroblast),
    npcs = 150
  )
erythroblast <- RunUMAP(erythroblast, dims = 1:100, min.dist = 0.4)
erythroblast <- FindNeighbors(erythroblast, dims = 1:100)
erythroblast <- FindClusters(erythroblast, resolution = 1)

erythroblast <- subset(erythroblast, subset = seurat_clusters != 1)
erythroblast <- subset(
  erythroblast,
  subset = (seurat_clusters == 3 &
    celltype == "hematopoietic precursor cell"),
  invert = TRUE
)
erythroblast <-
  RunPCA(erythroblast,
    features = VariableFeatures(object = erythroblast),
    npcs = 150
  )
erythroblast <- RunUMAP(erythroblast, dims = 1:100, min.dist = 0.3)
erythroblast <- FindNeighbors(erythroblast, dims = 1:100)
erythroblast <- FindClusters(erythroblast, resolution = 1)

# Monocle TI Monocytes
library(SeuratWrappers)
monocyte.cds <- as.cell_data_set(monocyte)
rownames(monocyte.cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
colnames(monocyte.cds@int_colData@listData$reducedDims@listData$UMAP) <-
  NULL
monocyte.cds <- learn_graph(monocyte.cds, use_partition = FALSE)
plot_cells(monocyte.cds,
  cell_size = 1,
  color_cells_by = "nCount_RNA"
)

plot_cells(
  monocyte.cds,
  color_cells_by = "celltype",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  cell_size = 1,
  group_label_size = 4,
  graph_label_size = 4
)

monocyte.cds <- order_cells(monocyte.cds)

monocyte <- AddMetaData(
  object = monocyte,
  metadata = monocyte.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

# Panel C
FeaturePlot(monocyte,
  features = c("Pseudotime"),
  pt.size = 0.05
) + scale_color_viridis_c() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Monocyte Differentiation")

# Panel D
FeaturePlot(monocyte,
  features = c("nCount_RNA"),
  pt.size = 0.05
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Transcripts")

pseudo_plot <- function(organ, gene, color, title) {
  organ_frame <- as.data.frame(organ@meta.data$Pseudotime)
  colnames(organ_frame) <- "Pseudotime"
  organ_frame[, gene] <- FetchData(organ, gene, cells = NULL, slot = "data")
  organ_frame[, "Transcripts"] <- FetchData(organ, "nCount_RNA", cells = NULL, slot = "data")
  organ_frame <-
    as.data.frame(organ_frame[is.finite(organ_frame$Pseudotime), ])
  colnames(organ_frame) <- c("Pseudotime", gene, "Transcripts")
  organ_frame$Pseudotime <- BBmisc::normalize(
    organ_frame$Pseudotime,
    method = "range",
    range = c(0, 1),
    margin = 2L
  )
  organ_frame$Transcripts <- BBmisc::normalize(
    organ_frame$Transcripts,
    method = "range",
    range = c(0, 2),
    margin = 2L
  )

  ggplot(organ_frame, aes(x = Pseudotime, y = eval(parse(text = (
    paste("`", gene, "`", sep = "")
  ))))) +
    theme_classic() +
    geom_point(color = "Grey", size = 1) +
    geom_smooth(
      method = "gam",
      color = color,
      size = 1.5,
      alpha = 1
    ) +
    theme(
      aspect.ratio = 0.6,
      text = element_text(size = 12),
      axis.text = element_text(size = 10, colour = "black"),
      plot.margin = margin(
        t = 7,
        r = 0,
        b = 7,
        l = 0,
        unit = "pt"
      ),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 11, hjust = 0.5)
    ) +
    labs(title = title, y = "Expression", x = "Pseudotime") +
    lims(x = c(0, 1)) +
    NoLegend()
}

# Panels E-H
pseudo_plot(monocyte, "nCount_RNA", "Black", "Transcripts")
pseudo_plot(monocyte, "Gapdh", "Black", "Gapdh")
pseudo_plot(monocyte, "Actb", "Black", "Actb")

pseudo_plot(monocyte, "Chd1", "maroon3", "Chd1")
pseudo_plot(monocyte, "Hira", "maroon3", "Hira")
pseudo_plot(monocyte, "Ino80", "maroon3", "Ino80")

pseudo_plot(monocyte, "Rpl9", "turquoise3", "Rpl9")
pseudo_plot(monocyte, "Rpl6", "turquoise3", "Rpl6")
pseudo_plot(monocyte, "Rps3", "turquoise3", "Rps3")

pseudo_plot(monocyte, "Ahnak", "chartreuse4", "Ahnak")
pseudo_plot(monocyte, "Mpeg1", "chartreuse4", "Mpeg1")
pseudo_plot(monocyte, "Emr1", "chartreuse4", "Emr1")

# Monocle TI Erythroblasts
library(SeuratWrappers)
erythroblast.cds <- as.cell_data_set(erythroblast)
rownames(erythroblast.cds@principal_graph_aux[["UMAP"]]$dp_mst) <-
  NULL
colnames(erythroblast.cds@int_colData@listData$reducedDims@listData$UMAP) <-
  NULL
erythroblast.cds <-
  learn_graph(erythroblast.cds, use_partition = FALSE)
plot_cells(erythroblast.cds,
  cell_size = 1,
  color_cells_by = "nCount_RNA"
)

plot_cells(
  erythroblast.cds,
  color_cells_by = "celltype",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  cell_size = 1,
  group_label_size = 4,
  graph_label_size = 4
)

erythroblast.cds <- order_cells(erythroblast.cds)

erythroblast <- AddMetaData(
  object = erythroblast,
  metadata = erythroblast.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

# Panel I
FeaturePlot(erythroblast,
  features = c("Pseudotime"),
  pt.size = 0.05
) + scale_color_viridis_c() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Erythroblast Differentiation") # 650 x 300

# Panel J
FeaturePlot(erythroblast,
  features = c("nCount_RNA"),
  pt.size = 0.05
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Transcripts") # 650 x 300

pseudo_plot <- function(organ, gene, color, title) {
  organ_frame <- as.data.frame(organ@meta.data$Pseudotime)
  colnames(organ_frame) <- "Pseudotime"
  organ_frame[, gene] <- FetchData(organ, gene, cells = NULL, slot = "data")
  organ_frame <-
    as.data.frame(organ_frame[is.finite(organ_frame$Pseudotime), ])
  colnames(organ_frame) <- c("Pseudotime", gene)
  organ_frame$Pseudotime <- BBmisc::normalize(
    organ_frame$Pseudotime,
    method = "range",
    range = c(0, 1),
    margin = 2L
  )

  ggplot(organ_frame, aes(x = Pseudotime, y = eval(parse(text = (
    paste("`", gene, "`", sep = "")
  ))))) +
    theme_classic() +
    geom_point(color = "Grey", size = 1) +
    geom_smooth(
      method = "gam",
      color = color,
      size = 1.5,
      alpha = 1
    ) +
    theme(
      aspect.ratio = 0.6,
      text = element_text(size = 12),
      axis.text = element_text(size = 10, colour = "black"),
      plot.margin = margin(
        t = 7,
        r = 0,
        b = 7,
        l = 0,
        unit = "pt"
      ),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 11, hjust = 0.5)
    ) +
    labs(title = title, y = "Expression", x = "Pseudotime") +
    lims(x = c(0, 1)) +
    NoLegend()
}
# Panels K-N
pseudo_plot(erythroblast, "nCount_RNA", "Black", "Transcripts")
pseudo_plot(erythroblast, "Gapdh", "Black", "Gapdh")
pseudo_plot(erythroblast, "Actb", "Black", "Actb")

pseudo_plot(erythroblast, "Chd1", "maroon3", "Chd1")
pseudo_plot(erythroblast, "Hira", "maroon3", "Hira")
pseudo_plot(erythroblast, "Ino80", "maroon3", "Ino80")

pseudo_plot(erythroblast, "Rpl9", "turquoise3", "Rpl9")
pseudo_plot(erythroblast, "Rpl6", "turquoise3", "Rpl6")
pseudo_plot(erythroblast, "Rps3", "turquoise3", "Rps3")

pseudo_plot(erythroblast, "Bpgm", "chartreuse4", "Bpgm")
pseudo_plot(erythroblast, "Beta-s", "chartreuse4", "Beta-s")
pseudo_plot(erythroblast, "Hbb-b2", "chartreuse4", "Hbb-b2")

# Panel O
erythroblast@assays$RNA@data <- as.matrix(log2(erythroblast@assays$RNA@counts +
  1))
erythroblast2 <- NormalizeData(erythroblast, normalization.method = "RC")
erythroblast2@assays$RNA@data <- as.matrix(log2(erythroblast2@assays$RNA@data +
  1))

pseudo_plot(erythroblast2, "Actb", "Black", "Actb")
pseudo_plot(erythroblast2, "Rps6", "turquoise3", "Rps6")
pseudo_plot(erythroblast2, "Chd1", "maroon3", "Chd1")
pseudo_plot(erythroblast2, "Beta-s", "chartreuse4", "Beta-s")
