library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)

tabulafacs <- readRDS("Tabula FACS/tabulafacsabs.rds")
tabulafacs@meta.data$logRNA <- log2(tabulafacs@meta.data$nCount_RNA + 1)
tabula_unnorm <- readRDS("Tabula FACS/tabulaunnorm.rds")
tabula_unnorm@meta.data$logRNA <- tabulafacs@meta.data$logRNA
tabula_unnorm@meta.data$facsRNA <- tabulafacs@meta.data$nCount_RNA
rm(tabulafacs)
cols <- read.csv(file = "Tabula FACS/tabulafacscols.csv")
rownames(cols) <- cols$Organs

autoprocess <- function(organ, tabula) {
  organ.seurat <- subset(tabula, subset = tissue == organ)
  organ.seurat <- NormalizeData(organ.seurat, normalization.method = "RC")
  organ.seurat@assays$RNA@data <- as.matrix(log2(organ.seurat@assays$RNA@counts +
    1))
  organ.seurat <-
    FindVariableFeatures(organ.seurat,
      selection.method = "vst",
      nfeatures = 2000
    )
  all.genes <- rownames(organ.seurat)
  organ.seurat <- ScaleData(organ.seurat, features = all.genes)
  organ.seurat <-
    RunPCA(organ.seurat, features = VariableFeatures(object = organ.seurat))
  organ.seurat <- FindNeighbors(organ.seurat, dims = 1:15)
  organ.seurat <- FindClusters(organ.seurat, resolution = 0.4)
  organ.seurat <-
    RunUMAP(organ.seurat, dims = 1:15, min.dist = 0.75)
  return(organ.seurat)
}

# auto_vector
auto_vector_log <- function(organ, tabula) {
  organ.seurat <- subset(tabula, subset = tissue == organ)
  organ.seurat <- NormalizeData(organ.seurat, normalization.method = "RC")
  organ.seurat@assays$RNA@data <- as.matrix(log2(organ.seurat@assays$RNA@counts +
    1))
  organ.seurat <-
    FindVariableFeatures(organ.seurat,
      selection.method = "vst",
      nfeatures = 5000
    )
  all.genes <- rownames(organ.seurat)
  organ.seurat <- ScaleData(organ.seurat, features = all.genes)
  organ.seurat <-
    RunPCA(organ.seurat, features = all.genes, npcs = 150)
  organ.seurat <-
    RunUMAP(organ.seurat, dims = 1:50, min.dist = 0.75)

  source("https://raw.githubusercontent.com/jumphone/Vector/master/Vector.R")
  VEC <- organ.seurat@reductions$umap@cell.embeddings
  rownames(VEC) <- colnames(organ.seurat)
  PCA <- organ.seurat@reductions$pca@cell.embeddings
  PCA <- vector.rankPCA(PCA)

  OUT <- vector.buildGrid(VEC, N = 150, SHOW = TRUE)
  OUT <- vector.buildNet(OUT, CUT = 1, SHOW = TRUE)
  OUT <- vector.getValue(OUT, PCA, SHOW = TRUE)
  organ.seurat@meta.data$QP <- OUT$VALUE
  return(organ.seurat)
}

# auto_QP_plot
auto_QP_plot <- function(organ_seurat, linecol, title) {
  organ_QP <- as.data.frame(organ_seurat@meta.data$logRNA)
  colnames(organ_QP) <- "RNA_Counts"
  organ_QP[, "QP"] <- organ_seurat@meta.data$QP

  ggplot(organ_QP, aes(x = QP, y = RNA_Counts)) +
    theme_classic() +
    geom_point(color = linecol, size = 1) +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(
        t = 0.5,
        r = 0.5,
        b = 0.5,
        l = 0.5,
        unit = "pt"
      )
    ) +
    labs(title = title, y = "", x = "") +
    stat_cor(
      label.x.npc = 0,
      label.y.npc = 1,
      method = "spearman",
      size = 4,
      fontface = "bold"
    ) +
    NoLegend()
}

# auto_feature_plot
auto_feature_plot <- function(organ_seurat, linecol, title) {
  organ_feature <- as.data.frame(organ_seurat@meta.data$logRNA)
  colnames(organ_feature) <- "RNA_Counts"
  organ_feature[, "feature"] <- as.data.frame(organ_seurat@meta.data$nFeature_Counts)

  ggplot(organ_feature, aes(x = feature, y = RNA_Counts)) +
    theme_classic() +
    geom_point(color = linecol, size = 1) +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(
        t = 0.5,
        r = 0.5,
        b = 0.5,
        l = 0.5,
        unit = "pt"
      )
    ) +
    labs(title = title, y = "", x = "") +
    stat_cor(
      label.x.npc = 0,
      label.y.npc = 1,
      method = "spearman",
      size = 4,
      fontface = "bold"
    ) +
    NoLegend() # 600 x 350
}

organlist <- cols$Organs
plotnumber <- 1
for (organ in organlist) {
  title <- gsub("_", " ", organ)
  organ_seurat <- auto_vector_log(organ, tabula_unnorm)
  eval(parse(
    text = paste(
      "pp",
      plotnumber,
      " = auto_feature_plot(organ_seurat, cols[organ, 'colors'], title)",
      sep = ""
    )
  ))
  plotnumber <- plotnumber + 1
  print(organ)
}

# Panel A, 1300 x 1100
library(egg)
gridpanels <- paste0("p", as.list(1:(plotnumber - 1)), collapse = ",")
eval(parse(text = paste(
  "grid.arrange(", gridpanels, ",ncol = 5)",
  sep = ""
)))

# Panel B, 1300 x 1100
gridpanels <- paste0("pp", as.list(1:(plotnumber - 1)), collapse = ",")
eval(parse(text = paste(
  "grid.arrange(", gridpanels, ",ncol = 5)",
  sep = ""
)))
