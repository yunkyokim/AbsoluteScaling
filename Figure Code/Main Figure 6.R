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
granulocyte <- subset(marrow, subset = celltype == "hematopoietic precursor cell" |
  celltype == "granulocytopoietic cell" |
  celltype == "granulocyte")

granulocyte@assays$RNA@data <- as.matrix(log2(granulocyte@assays$RNA@counts + 1))
granulocyte <-
  FindVariableFeatures(granulocyte,
    selection.method = "vst",
    nfeatures = 5000
  )

all.genes <- rownames(granulocyte)
granulocyte <- ScaleData(granulocyte, features = all.genes)
granulocyte <-
  RunPCA(granulocyte,
    features = VariableFeatures(object = granulocyte),
    npcs = 150
  )

granulocyte <- RunUMAP(granulocyte, dims = 1:60, min.dist = 0.4)
granulocyte <- FindNeighbors(granulocyte, dims = 1:60)
granulocyte <- FindClusters(granulocyte)

granulocyte <- subset(granulocyte, subset = seurat_clusters != 6 &
  seurat_clusters != 9)
granulocyte <- RunUMAP(granulocyte, dims = 1:60, min.dist = 0.4)
granulocyte <- FindNeighbors(granulocyte, dims = 1:60)
granulocyte <- FindClusters(granulocyte)


# Monocle TI
library(SeuratWrappers)
granulocyte.cds <- as.cell_data_set(granulocyte)
granulocyte.cds <-
  cluster_cells(cds = granulocyte.cds, reduction_method = "UMAP")
granulocyte.cds <-
  learn_graph(granulocyte.cds, use_partition = TRUE)
plot_cells(granulocyte.cds,
  cell_size = 1,
  color_cells_by = "nCount_RNA"
)

plot_cells(
  granulocyte.cds,
  color_cells_by = "celltype",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  cell_size = 1,
  group_label_size = 4,
  graph_label_size = 4
)

rownames(granulocyte.cds@principal_graph_aux[["UMAP"]]$dp_mst) <-
  NULL
colnames(granulocyte.cds@int_colData@listData$reducedDims@listData$UMAP) <-
  NULL
granulocyte.cds <- order_cells(granulocyte.cds)

granulocyte <- AddMetaData(
  object = granulocyte,
  metadata = granulocyte.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

# Panel A
FeaturePlot(granulocyte,
  features = c("Pseudotime"),
  pt.size = 0.05
) + scale_color_viridis_c() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Granulocyte Differentiation") # 650 x 350

# Panel B
FeaturePlot(granulocyte,
  features = c("nCount_RNA"),
  pt.size = 0.05
) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "Transcripts")

Idents(granulocyte) <- "celltype"
granulocyte_markers <- FindMarkers(granulocyte, ident.1 = "hematopoietic precursor cell", ident.2 = "granulocyte")

mat <- as.data.frame(granulocyte@assays$RNA@data)
mat["pseudotime", ] <- (granulocyte@meta.data$Pseudotime)
mat <- as.data.frame(t(mat))
mat <- mat[order(mat$pseudotime), ]
mat$pseudotime <- NULL
mat <-
  mat[, (which(colnames(mat) %in% rownames(granulocyte_markers)))]

mat2 <- mat
mat2 <- BBmisc::normalize(mat2,
  method = "range",
  range = c(0, 2),
  margin = 2L
)
mat2 <- as.matrix(t(mat2))

library(RColorBrewer)
newcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))
puor <- newcol(200)

library(ComplexHeatmap)
HM <- Heatmap(
  mat2,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  col = puor,
  column_title = "Pseudotime",
  show_row_dend = FALSE,
  km = 6,
  name = "Expression",
  use_raster = TRUE,
  right_annotation = ha
)

# Panel C, 700 x 400
HM <- draw(HM)


cluster_list <- list(
  text1 = "Class 1",
  text2 = "Class 2",
  text3 = "Class 3",
  text4 = "Class 4",
  text5 = "Class 5",
  text5 = "Class 6"
)

ha <- rowAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(unlist(cluster_list)) + unit(4, "mm")))

for (i in 1:6) {
  decorate_annotation("foo", slice = i, {
    grid.rect(
      x = 0,
      width = unit(2, "mm"),
      gp = gpar(fill = i + 1, col = NA),
      just = "left"
    )
    grid.text(paste(cluster_list[[i]], collapse = "\n"),
      x = unit(4, "mm"),
      just = "left"
    )
  })
}

# Cluster Extraction
rcl.list <- row_order(HM)
lapply(rcl.list, function(x) {
  length(x)
})

library(magrittr)
clu_df <- lapply(names(rcl.list), function(i) {
  out <-
    data.frame(
      GeneID = rownames(mat2[rcl.list[[i]], ]),
      Cluster = paste0("cluster", i),
      stringsAsFactors = FALSE
    )
  return(out)
}) %>%
  do.call(rbind, .)
write.csv(clu_df, "Tabula 10X BM/granulocyte_pseudotime_clusters.csv")

mat <- as.data.frame(granulocyte@assays$RNA@data)
mat <- as.data.frame(t(mat))
mat <-
  mat[, (which(colnames(mat) %in% rownames(granulocyte_markers)))]

mat3 <- as.data.frame(t(mat))
mat3$gene <- rownames(mat3)
mat3[, "clusters"] <- clu_df$Cluster[match(mat3$gene, clu_df$GeneID)]
mat3$gene <- NULL
mat3 <- aggregate(. ~ clusters, mat3, mean)
rownames(mat3) <- mat3$clusters
mat3$clusters <- NULL
mat3 <- as.data.frame(t(mat3))
mat3 <- BBmisc::normalize(mat3,
  method = "range",
  range = c(0, 2),
  margin = 2L
)

granulocyte@meta.data$cluster1 <- mat3$cluster1
granulocyte@meta.data$cluster2 <- mat3$cluster2
granulocyte@meta.data$cluster3 <- mat3$cluster3
granulocyte@meta.data$cluster4 <- mat3$cluster4
granulocyte@meta.data$cluster5 <- mat3$cluster5
granulocyte@meta.data$cluster6 <- mat3$cluster6

mat3$average <- granulocyte$nCount_RNA
mat3$average <- log2(mat3$average + 1)
mat3$average <- BBmisc::normalize(
  mat3$average,
  method = "range",
  range = c(0, 2),
  margin = 2L
)

# Plotting
pseudo_plot <- function(organ, gene, color, title, fill) {
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
    paste(gene, sep = "")
  ))))) +
    theme_classic() +
    geom_point(color = "Grey", size = 1) +
    geom_smooth(
      method = "loess",
      color = color,
      size = 1.5,
      alpha = 1,
      se = TRUE,
      level = 0.99,
      fill = fill
    ) +
    geom_smooth(
      aes(y = Transcripts),
      color = "black",
      linetype = "twodash",
      se = TRUE,
      level = 0.99
    ) +
    theme(
      aspect.ratio = 1.2,
      text = element_text(size = 12),
      axis.text = element_text(size = 10, colour = "black"),
      plot.margin = margin(
        t = 7,
        r = 0,
        b = 7,
        l = 0,
        unit = "pt"
      ),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 11, hjust = 0.5)
    ) +
    labs(title = title, y = "Expression", x = "Pseudotime") +
    lims(x = c(0, 1)) +
    NoLegend()
}

# Panels E-I Left
pseudo_plot(
  granulocyte,
  "cluster1",
  "#DF536B",
  "Class 1 Transcripts",
  "#eb97a6"
)
pseudo_plot(
  granulocyte,
  "cluster2",
  "#61D04F",
  "Class 2 Transcripts",
  "#a0e295"
)
pseudo_plot(
  granulocyte,
  "cluster3",
  "#2297E6",
  "Class 3 Transcripts",
  "#7ac0f0"
)
pseudo_plot(
  granulocyte,
  "cluster4",
  "#28E2E5",
  "Class 4 Transcripts",
  "#7eedef"
)
pseudo_plot(
  granulocyte,
  "cluster5",
  "#CD0BBC",
  "Class 5 Transcripts",
  "#e16cd6"
)
pseudo_plot(
  granulocyte,
  "cluster6",
  "#F5C710",
  "Class 6 Transcripts",
  "#f9dd6f"
)

# DAVID analysis done externally based on cluster marker analysis
auto_goplot <- function(filepath, color) {
  david <- read.delim(filepath)
  david$logp <- -log(david$PValue)
  david$Name <- sapply(X = strsplit(david$Term, split = "~"), FUN = "[", 2)

  print(
    ggplot(head(david, 10), aes(x = logp, y = reorder(Name, logp)), ) +
      geom_bar(stat = "identity", fill = color) +
      theme_classic() +
      NoLegend() +
      theme(
        plot.margin = margin(
          t = 5,
          r = 5,
          b = 5,
          l = 0,
          unit = "pt"
        ),
        text = element_text(size = 11),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1.5,
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
      ) +
      labs(title = "GO BP Terms", y = "", x = "-log p-Value") +
      scale_y_discrete(
        label = function(x) {
          stringr::str_trunc(x, 30)
        }
      )
  )
}

# Panel D-I Right
auto_goplot("granulocytecluster1.txt", "#DF536B")
auto_goplot("granulocytecluster2.txt", "#61D04F")
auto_goplot("granulocytecluster3.txt", "#2297E6")
auto_goplot("granulocytecluster4.txt", "#28E2E5")
auto_goplot("granulocytecluster5.txt", "#CD0BBC")
auto_goplot("granulocytecluster6.txt", "#F5C710")
