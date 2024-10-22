library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(scran)
library(scRNAseq)
library(scater)
library(DoubletFinder)

# setwd("current_dir/")
muscle <- readRDS("Oprescu2020/muscle.rds")

doublet_detection <- function(seurat_dataset, est_doublet) {
  # Log-Normalized Doublet Detection
  seurat_dataset <- NormalizeData(seurat_dataset, normalization.method = "LogNormalize")
  seurat_dataset <- ScaleData(seurat_dataset)
  seurat_dataset <-
    FindVariableFeatures(seurat_dataset,
      selection.method = "vst",
      nfeatures = 2000
    )
  seurat_dataset <- RunPCA(seurat_dataset)
  seurat_dataset <- FindNeighbors(seurat_dataset, dims = 1:10)
  seurat_dataset <- FindClusters(seurat_dataset, resolution = 0.4)
  seurat_dataset <- RunUMAP(seurat_dataset, dims = 1:10)

  # Doublet Detection
  sweep.res <- paramSweep_v3(seurat_dataset, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  log.bcmvn <- (as.data.frame(find.pK(sweep.stats)))
  log.bcmvn$pK <- as.numeric(as.character(log.bcmvn$pK))
  max_pk <- log.bcmvn$pK[[which(grepl(max(data.matrix(
    log.bcmvn$BCmetric
  )), log.bcmvn$BCmetric))]]

  annotations <- seurat_dataset@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  est_doublet <- 0.04
  nExp_poi <- round(est_doublet * length(seurat_dataset$orig.ident))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  seurat_dataset <-
    doubletFinder_v3(
      seurat_dataset,
      PCs = 1:10,
      pN = 0.25,
      pK = max_pk,
      nExp = nExp_poi.adj,
      reuse.pANN = FALSE,
      sct = FALSE
    )
  meta_doublets <- paste(
    "DF.classifications_0.25_",
    as.character(max_pk),
    "_",
    as.character(nExp_poi.adj),
    sep = ""
  )

  # Subset
  seurat_dataset <- eval(parse(
    text = paste(
      "subset(seurat_dataset, subset = ",
      meta_doublets,
      "== 'Singlet')",
      sep = ""
    )
  ))
  rm(
    log.bcmvn,
    sweep.res,
    sweep.stats,
    annotations,
    homotypic.prop,
    nExp_poi,
    nExp_poi.adj,
    max_pk,
    est_doublet,
    meta_doublets
  )
  return(seurat_dataset)
}
muscleregen <- doublet_detection(muscle, 0.04)
muscleregen@assays$RNA@data <- as.matrix(log2(muscleregen@assays$RNA@counts +
  1))
muscleregen <- FindVariableFeatures(muscleregen)
muscleregen <- ScaleData(muscleregen)
muscleregen <-
  RunPCA(
    muscleregen,
    npcs = 100,
    ndims.print = 1:5,
    nfeatures.print = 5
  )
muscleregen <- FindNeighbors(muscleregen, dims = 1:15)
muscleregen <- FindClusters(muscleregen, resolution = 2)
muscleregen <- RunUMAP(muscleregen, dims = 1:15, min.dist = 0.75)

# Figure Panels
Idents(muscleregen) <- "metacluster"
metanames <- muscleregen@meta.data$metacluster
metanames <- gsub("_", " ", metanames)
muscleregen@meta.data$metacluster <- metanames
highlight <- WhichCells(muscleregen, idents = c("MuSC"))

# Panel A, 340 x 350
DimPlot(
  muscleregen,
  label = F,
  cells.highlight = highlight,
  sizes.highlight = 0.1,
  cols.highlight = c("turquoise3"),
  cols = "grey",
  pt.size = 0.1
) +
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
  ) + NoLegend()

# Panel C, 600 x 350
vlncols <- c("turquoise3", rep("grey", 14))
VlnPlot(
  muscleregen,
  features = c("nCount_RNA"),
  cols = vlncols,
  ncol = 1,
  pt.size = 0,
  sort = TRUE
) +
  ylab("Absolute UMI") + xlab("") + NoLegend() +
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 8,
    colour = "grey27",
    shape = 95
  ) +
  scale_y_continuous(
    labels = function(x) {
      format(x, scientific = TRUE)
    }
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
    )
  )

# MuSC Cell Subset
musc <- subset(muscleregen, subset = metacluster == "MuSC")
musc@assays$RNA@data <- as.matrix(log2(musc@assays$RNA@counts + 1))
musc <- FindVariableFeatures(musc)
musc <- ScaleData(musc)
musc <-
  RunPCA(
    musc,
    npcs = 100,
    ndims.print = 1:5,
    nfeatures.print = 5
  )
musc <- FindNeighbors(musc, dims = 1:15)
musc <- FindClusters(musc, resolution = 2)
musc <- RunUMAP(musc, dims = 1:15, min.dist = 0.75)

Idents(musc) <- "timepoint"
mlevels <-
  c(
    "Noninjured",
    "0.5 DPI",
    "2 DPI",
    "3.5 DPI",
    "5 DPI",
    "10 DPI",
    "21 DPI"
  )
musc@active.ident <-
  factor(
    x = musc@active.ident,
    levels = c(
      "Noninjured",
      "0.5 DPI",
      "2 DPI",
      "3.5 DPI",
      "5 DPI",
      "10 DPI",
      "21 DPI"
    )
  )
colfunc <- colorRampPalette(c("turquoise4", "turquoise2"))

# Panel B, 600 x 350
DimPlot(musc, pt.size = 1) +
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

musc@meta.data$logRNA <- log2(musc@meta.data$nCount_RNA + 1)
cellstates <- as.data.frame(table(as.data.frame(musc@active.ident)))
celltranscripts <- data.frame(
  Transcripts = musc@meta.data$nCount_RNA,
  timepoint = musc@meta.data$timepoint
)
celltranscripts <- aggregate(. ~ timepoint, celltranscripts, median)
cellstates$transcripts <- celltranscripts$Transcripts[match(cellstates$Var1, celltranscripts$timepoint)]
celltranscripts <- data.frame(
  Transcripts = musc@meta.data$nCount_RNA,
  timepoint = musc@meta.data$timepoint
)
celltranscripts <- aggregate(. ~ timepoint, celltranscripts, sd)
cellstates$stdev <- celltranscripts$Transcripts[match(cellstates$Var1, celltranscripts$timepoint)]
cellstates$colors <- colfunc(7)

# Panel D
ggplot(cellstates, aes(x = Var1)) +
  geom_line(aes(y = Freq, group = 1, colour = "Cell Number"), size = 2) +
  geom_line(aes(
    y = transcripts / 140,
    group = 1,
    colour = "Transcripts"
  ), size = 2) +
  geom_errorbar(
    aes(
      ymin = (transcripts - stdev) / 140,
      ymax = (transcripts + stdev) / 140
    ),
    width = .2,
    position = position_dodge(.9)
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 0.5,
    text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y.right = element_text(margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 10
    ))
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * 140, name = "Transcripts (UMIs)"),
    name = "Cell Number"
  ) +
  labs(title = "", y = "", x = "") +
  scale_colour_manual(values = c("turquoise4", "grey"), name = "")

entrymarkers <- FindMarkers(musc, ident.1 = "0.5 DPI", ident.2 = "Noninjured")
write.table(entrymarkers, "markers.csv", sep = ",")
goterms <- read.csv("Oprescu2020/goterms.csv")
goterms$logp <- -log(goterms$Adjusted.P.value)
goterms$Name <- sapply(
  X = strsplit(goterms$ï..Term, split = " \\(GO"),
  FUN = "[",
  1
)

# Panel E, 1000 x 830
ggplot(head(goterms, 15), aes(x = logp, y = reorder(Name, logp)), ) +
  geom_bar(stat = "identity", fill = "turquoise4") +
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
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 1),
    aspect.ratio = 1,
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  xlim(0, 45) +
  labs(title = "", y = "", x = "-log p-Value") +
  scale_y_discrete(
    label = function(x) {
      stringr::str_trunc(x, 60)
    }
  )
