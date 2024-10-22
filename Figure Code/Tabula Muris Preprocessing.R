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

# Tabula FACS Preprocessing
tabula_raw <- readRDS("Tabula FACS/tabularaw.rds") # Merged object from raw read files of original Tabula Muris

# ERCC Normalized
filteredcounts <- GetAssayData(object = tabulafacs, assay = "RNA", slot = "counts")
scran.data <- SingleCellExperiment(assays = list(counts = filteredcounts))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <- splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))
altExpNames(scran.data)
scran.data <- computeSpikeFactors(scran.data, "ERCC")
summary(sizeFactors(scran.data))
scran.data <- logNormCounts(scran.data, log = FALSE)
tabulafacs <- CreateSeuratObject(counts = as.sparse(x = assay(scran.data, "normcounts")))

tabulafacs <- AddMetaData(tabulafacs, PercentageFeatureSet(tabulafacs, pattern = "^ERCC-"), col.name = "ERCC")
tabulafacs@meta.data$noclusters <- 0
tabulafacs <- subset(tabulafacs, subset = nCount_RNA < 7500000)

# Annotation
tabulamat <- as.data.frame(tabulafacs@assays$RNA@counts)
cellbarcodes <- sapply(X = strsplit(colnames(tabulamat), split = "_"), FUN = "[", 2)
colnames(tabulamat) <- cellbarcodes
tabulamat <- tabulamat[, order(colnames(tabulamat))]
head(tabulamat)

annotations <- as.data.frame(read.csv("Tabula FACS/annotations_facs.csv"), sep = ",", header = TRUE)
annotations$cell <- sapply(X = strsplit(annotations$cell, split = "_"), FUN = "[", 1)
annotations <- annotations[annotations$cell %in% cellbarcodes, ]
annotations <- annotations[order(annotations$cell), ]

tabulacol <- colnames(tabulamat)
tabulamat <- tabulamat[, tabulacol %in% annotations$cell]
tabulamat <- tabulamat[, order(colnames(tabulamat))]
annotations <- annotations[order(annotations$cell), ]
rm(tabulacol)

tabulafacs <- CreateSeuratObject(counts = tabulamat, project = "tabulafacs")
tabulafacs@meta.data$tissue <- annotations$tissue
tabulafacs@meta.data$mouse_sex <- annotations$mouse.sex
tabulafacs@meta.data$celltype <- annotations$cell_ontology_class
tabulafacs@meta.data$celltypeid <- annotations$cell_ontology_id
tabulafacs@meta.data$noclusters <- 0

saveRDS(tabulafacs, "Tabula FACS/tabulafacsabs.rds")


# Global scaling
tabula_log <- tabula_raw

filteredcounts <- GetAssayData(object = tabula_log, assay = "RNA", slot = "counts")
scran.data <- SingleCellExperiment(assays = list(counts = filteredcounts))
is.spike <- grepl("^ERCC-", rownames(scran.data))
scran.data <- splitAltExps(scran.data, ifelse(is.spike, "ERCC", "gene"))
altExpNames(scran.data)
tabula_log <- CreateSeuratObject(counts = as.sparse(x = assay(scran.data, "counts")))

tabula_log <- AddMetaData(tabula_log, PercentageFeatureSet(tabula_log, pattern = "^ERCC-"), col.name = "ERCC")
tabula_log@meta.data$noclusters <- 0

# Annotation
tabulamat <- as.data.frame(tabula_log@assays$RNA@counts)
cellbarcodes <- sapply(X = strsplit(colnames(tabulamat), split = "_"), FUN = "[", 2)
colnames(tabulamat) <- cellbarcodes
tabulamat <- tabulamat[, order(colnames(tabulamat))]

annotations <- as.data.frame(read.csv("Tabula FACS/annotations_facs.csv"), sep = ",", header = TRUE)
annotations$cell <- sapply(X = strsplit(annotations$cell, split = "_"), FUN = "[", 1)
annotations <- annotations[annotations$cell %in% cellbarcodes, ]
annotations <- annotations[order(annotations$cell), ]

tabulacol <- colnames(tabulamat)
tabulamat <- tabulamat[, tabulacol %in% annotations$cell]

tabulaercclist <- colnames(tabulafacs@assays$RNA@counts)

tabulamat <- tabulamat[, intersect(tabulacol, tabulaercclist)]
annotations <- annotations[intersect(annotations$cell, tabulaercclist), ]

tabulamat <- tabulamat[, order(colnames(tabulamat))]
annotations <- annotations[order(annotations$cell), ]
rm(tabulacol)

tabula_log <- CreateSeuratObject(counts = tabulamat, project = "tabula_log")
tabula_log@meta.data$tissue <- annotations$tissue
tabula_log@meta.data$mouse_sex <- annotations$mouse.sex
tabula_log@meta.data$celltype <- annotations$cell_ontology_class
tabula_log@meta.data$celltypeid <- annotations$cell_ontology_id
tabula_log@meta.data$noclusters <- 0

saveRDS(tabula_log, "Tabula FACS/tabulafacsglo.rds")


# No Scaling
tabula_unnorm <- tabula_raw

# Annotation
tabulamat <- as.data.frame(tabula_unnorm@assays$RNA@counts)
cellbarcodes <- sapply(X = strsplit(colnames(tabulamat), split = "_"), FUN = "[", 2)
colnames(tabulamat) <- cellbarcodes
tabulamat <- tabulamat[, order(colnames(tabulamat))]

annotations <- as.data.frame(read.csv("Tabula FACS/annotations_facs.csv"), sep = ",", header = TRUE)
annotations$cell <- sapply(X = strsplit(annotations$cell, split = "_"), FUN = "[", 1)
annotations <- annotations[annotations$cell %in% cellbarcodes, ]
annotations <- annotations[order(annotations$cell), ]

tabulacol <- colnames(tabulamat)
tabulamat <- tabulamat[, tabulacol %in% annotations$cell]

tabulaercclist <- colnames(tabulafacs@assays$RNA@counts)

tabulamat <- tabulamat[, intersect(tabulacol, tabulaercclist)]
rownames(annotations) <- annotations$cell
annotations <- annotations[intersect(annotations$cell, tabulaercclist), ]

tabulamat <- tabulamat[, order(colnames(tabulamat))]
annotations <- annotations[order(annotations$cell), ]
rm(tabulacol)

tabula_unnorm <- CreateSeuratObject(counts = tabulamat, project = "tabula_unnorm")
tabula_unnorm@meta.data$tissue <- annotations$tissue
tabula_unnorm@meta.data$mouse_sex <- annotations$mouse.sex
tabula_unnorm@meta.data$celltype <- annotations$cell_ontology_class
tabula_unnorm@meta.data$celltypeid <- annotations$cell_ontology_id
tabula_unnorm@meta.data$noclusters <- 0
tabula_unnorm <- AddMetaData(tabula_unnorm, PercentageFeatureSet(tabula_unnorm, pattern = "^ERCC-"), col.name = "ERCC")

saveRDS(tabula_unnorm, "Tabula FACS/tabulaunnorm.rds")
rm(list = setdiff(ls(), c("tabulafacs", "tabula_log", "tabula_raw", "tabula_unnorm")))

# Tabula 10X
# Import and Processing
bladder1 <- Read10X(data.dir = "Tabula 10X/Bladder-10X_P4_3/")
bladder2 <- Read10X(data.dir = "Tabula 10X/Bladder-10X_P4_4/")
bladder3 <- Read10X(data.dir = "Tabula 10X/Bladder-10X_P7_7/")
bladder1.data <- CreateSeuratObject(counts = bladder1, project = "bladder1", min.cells = 3, min.features = 200)
bladder2.data <- CreateSeuratObject(counts = bladder2, project = "bladder2", min.cells = 3, min.features = 200)
bladder3.data <- CreateSeuratObject(counts = bladder3, project = "bladder3", min.cells = 3, min.features = 200)
bladder.data <- merge(bladder1.data, y = c(bladder2.data, bladder3.data), add.cell.ids = c("10X_P4_3", "10X_P4_4", "10X_P7_7"), project = "bladder")

aorta <- Read10X(data.dir = "Tabula 10X/Heart_and_Aorta-10X_P7_4/")
aorta.data <- CreateSeuratObject(counts = aorta, project = "aorta", min.cells = 3, min.features = 200)
aorta.data <- RenameCells(object = aorta.data, add.cell.id = "10X_P7_4")

kidney1 <- Read10X(data.dir = "Tabula 10X/Kidney-10X_P4_5/")
kidney2 <- Read10X(data.dir = "Tabula 10X/Kidney-10X_P4_6/")
kidney3 <- Read10X(data.dir = "Tabula 10X/Kidney-10X_P7_5/")
kidney1.data <- CreateSeuratObject(counts = kidney1, project = "kidney1", min.cells = 3, min.features = 200)
kidney2.data <- CreateSeuratObject(counts = kidney2, project = "kidney2", min.cells = 3, min.features = 200)
kidney3.data <- CreateSeuratObject(counts = kidney3, project = "kidney3", min.cells = 3, min.features = 200)
kidney.data <- merge(kidney1.data, y = c(kidney2.data, kidney3.data), add.cell.ids = c("10X_P4_5", "10X_P4_6", "10X_P7_5"), project = "kidney")

limb1 <- Read10X(data.dir = "Tabula 10X/Limb_Muscle-10X_P7_14/")
limb2 <- Read10X(data.dir = "Tabula 10X/Limb_Muscle-10X_P7_15/")
limb1.data <- CreateSeuratObject(counts = limb1, project = "limb1", min.cells = 3, min.features = 200)
limb2.data <- CreateSeuratObject(counts = limb2, project = "limb2", min.cells = 3, min.features = 200)
limb.data <- merge(limb1.data, y = c(limb2.data), add.cell.ids = c("10X_P7_14", "10X_P7_15"), project = "limb")

liver1 <- Read10X(data.dir = "Tabula 10X/Liver-10X_P4_2/")
liver2 <- Read10X(data.dir = "Tabula 10X/Liver-10X_P7_0/")
liver3 <- Read10X(data.dir = "Tabula 10X/Liver-10X_P7_1/")
liver1.data <- CreateSeuratObject(counts = liver1, project = "liver1", min.cells = 3, min.features = 200)
liver2.data <- CreateSeuratObject(counts = liver2, project = "liver2", min.cells = 3, min.features = 200)
liver3.data <- CreateSeuratObject(counts = liver3, project = "liver3", min.cells = 3, min.features = 200)
liver.data <- merge(liver1.data, y = c(liver2.data, liver3.data), add.cell.ids = c("10X_P4_2", "10X_P7_0", "10X_P7_1"), project = "liver")

lung1 <- Read10X(data.dir = "Tabula 10X/Lung-10X_P7_8/")
lung2 <- Read10X(data.dir = "Tabula 10X/Lung-10X_P7_9/")
lung3 <- Read10X(data.dir = "Tabula 10X/Lung-10X_P8_12/")
lung4 <- Read10X(data.dir = "Tabula 10X/Lung-10X_P8_13/")
lung1.data <- CreateSeuratObject(counts = lung1, project = "lung1", min.cells = 3, min.features = 200)
lung2.data <- CreateSeuratObject(counts = lung2, project = "lung2", min.cells = 3, min.features = 200)
lung3.data <- CreateSeuratObject(counts = lung3, project = "lung3", min.cells = 3, min.features = 200)
lung4.data <- CreateSeuratObject(counts = lung4, project = "lung4", min.cells = 3, min.features = 200)
lung.data <- merge(lung1.data, y = c(lung2.data, lung3.data, lung4.data), add.cell.ids = c("10X_P7_8", "10X_P7_9", "10X_P8_12", "10X_P8_13"), project = "lung")

mammary1 <- Read10X(data.dir = "Tabula 10X/Mammary_Gland-10X_P7_12/")
mammary2 <- Read10X(data.dir = "Tabula 10X/Mammary_Gland-10X_P7_13/")
mammary1.data <- CreateSeuratObject(counts = mammary1, project = "mammary1", min.cells = 3, min.features = 200)
mammary2.data <- CreateSeuratObject(counts = mammary2, project = "mammary2", min.cells = 3, min.features = 200)
mammary.data <- merge(mammary1.data, y = c(mammary2.data), add.cell.ids = c("10X_P7_12", "10X_P7_13"), project = "mammary")

marrow1 <- Read10X(data.dir = "Tabula 10X/Marrow-10X_P7_2/")
marrow2 <- Read10X(data.dir = "Tabula 10X/Marrow-10X_P7_3/")
marrow1.data <- CreateSeuratObject(counts = marrow1, project = "marrow1", min.cells = 3, min.features = 200)
marrow2.data <- CreateSeuratObject(counts = marrow2, project = "marrow2", min.cells = 3, min.features = 200)
marrow.data <- merge(marrow1.data, y = c(marrow2.data), add.cell.ids = c("10X_P7_2", "10X_P7_3"), project = "marrow")

spleen1 <- Read10X(data.dir = "Tabula 10X/Spleen-10X_P4_7/")
spleen2 <- Read10X(data.dir = "Tabula 10X/Spleen-10X_P7_6/")
spleen1.data <- CreateSeuratObject(counts = spleen1, project = "spleen1", min.cells = 3, min.features = 200)
spleen2.data <- CreateSeuratObject(counts = spleen2, project = "spleen2", min.cells = 3, min.features = 200)
spleen.data <- merge(spleen1.data, y = c(spleen2.data), add.cell.ids = c("10X_P4_7", "10X_P7_6"), project = "spleen")

thymus <- Read10X(data.dir = "Tabula 10X/Thymus-10X_P7_11/")
thymus.data <- CreateSeuratObject(counts = thymus, project = "thymus", min.cells = 3, min.features = 200)
thymus.data <- RenameCells(object = thymus.data, add.cell.id = "10X_P7_11")

tongue1 <- Read10X(data.dir = "Tabula 10X/Tongue-10X_P4_0/")
tongue2 <- Read10X(data.dir = "Tabula 10X/Tongue-10X_P4_1/")
tongue3 <- Read10X(data.dir = "Tabula 10X/Tongue-10X_P7_10/")
tongue1.data <- CreateSeuratObject(counts = tongue1, project = "tongue1", min.cells = 3, min.features = 200)
tongue2.data <- CreateSeuratObject(counts = tongue2, project = "tongue2", min.cells = 3, min.features = 200)
tongue3.data <- CreateSeuratObject(counts = tongue3, project = "tongue3", min.cells = 3, min.features = 200)
tongue.data <- merge(tongue1.data, y = c(tongue2.data, tongue3.data), add.cell.ids = c("10X_P4_0", "10X_P4_1", "10X_P7_10"), project = "tongue")

trachea1 <- Read10X(data.dir = "Tabula 10X/Trachea-10X_P8_14/")
trachea2 <- Read10X(data.dir = "Tabula 10X/Trachea-10X_P8_15/")
trachea1.data <- CreateSeuratObject(counts = trachea1, project = "trachea1", min.cells = 3, min.features = 200)
trachea2.data <- CreateSeuratObject(counts = trachea2, project = "trachea2", min.cells = 3, min.features = 200)
trachea.data <- merge(trachea1.data, y = c(trachea2.data), add.cell.ids = c("10X_P8_14", "10X_P8_15"), project = "trachea")

# Doublet Detection
rm(list = ls()[!(ls() %in% c(
  "aorta.data", "bladder.data", "kidney.data", "limb.data",
  "liver.data", "lung.data", "mammary.data", "marrow.data", "spleen.data", "thymus.data", "tongue.data", "trachea.data"
))])
organs.data <- as.list(ls())

auto_doublet_detection <- function(seurat_dataset, est_doublet) {
  # Log-Normalized Doublet Detection
  seurat_dataset <- subset(seurat_dataset, subset = nFeature_RNA > 1000 & nCount_RNA > 500 & nCount_RNA < 60000)
  seurat_dataset <- NormalizeData(seurat_dataset, normalization.method = "LogNormalize")
  seurat_dataset <- ScaleData(seurat_dataset)
  seurat_dataset <- FindVariableFeatures(seurat_dataset, selection.method = "vst", nfeatures = 2000)
  seurat_dataset <- RunPCA(seurat_dataset)
  seurat_dataset <- FindNeighbors(seurat_dataset, dims = 1:10)
  seurat_dataset <- FindClusters(seurat_dataset, resolution = 0.4)
  seurat_dataset <- RunUMAP(seurat_dataset, dims = 1:10)

  # Doublet Detection
  sweep.res <- paramSweep_v3(seurat_dataset, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  log.bcmvn <- (as.data.frame(find.pK(sweep.stats)))
  log.bcmvn$pK <- as.numeric(as.character(log.bcmvn$pK))
  max_pk <- log.bcmvn$pK[[which(grepl(max(data.matrix(log.bcmvn$BCmetric)), log.bcmvn$BCmetric))]]

  annotations <- seurat_dataset@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  est_doublet <- 0.04
  nExp_poi <- round(est_doublet * length(seurat_dataset$orig.ident))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  seurat_dataset <- doubletFinder_v3(seurat_dataset, PCs = 1:10, pN = 0.25, pK = max_pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  meta_doublets <- paste("DF.classifications_0.25_", as.character(max_pk), "_", as.character(nExp_poi.adj), sep = "")

  # Subset
  seurat_dataset <- eval(parse(text = paste("subset(seurat_dataset, subset = ", meta_doublets, "== 'Singlet')", sep = "")))
  rm(log.bcmvn, sweep.res, sweep.stats, annotations, homotypic.prop, nExp_poi, nExp_poi.adj, max_pk, est_doublet, meta_doublets)
  return(seurat_dataset)
}

for (organ_dataset in organs.data) {
  eval(parse(text = paste(organ_dataset, "= auto_doublet_detection(", organ_dataset, ", 0.04)", sep = "")))
}

# Merge
tabula10X <- merge(aorta.data, y = c(
  bladder.data,
  kidney.data,
  limb.data,
  liver.data,
  lung.data,
  mammary.data,
  marrow.data,
  spleen.data,
  thymus.data,
  tongue.data,
  trachea.data
), project = "tabula10X")

rm(list = setdiff(ls(), "tabula10X"))

tabula10X@meta.data$noclusters <- 0
tabula10X <- subset(tabula10X, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & nCount_RNA < 60000)

tabula10X <- FindVariableFeatures(tabula10X)
all.genes <- rownames(tabula10X)
tabula10X <- ScaleData(tabula10X, features = all.genes)
tabula10X <- RunPCA(tabula10X, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

tabula10X <- FindNeighbors(tabula10X, dims = 1:15)
tabula10X <- FindClusters(tabula10X, resolution = 2)
tabula10X <- RunUMAP(tabula10X, dims = 1:15, min.dist = 0.75)

# Annotation
tabulamat <- as.data.frame(tabula10X@assays$RNA@counts)
cellbarcodes <- sapply(X = strsplit(colnames(tabulamat), split = "-"), "[", 1) #####
colnames(tabulamat) <- cellbarcodes

annotations <- as.data.frame(read.csv(file = paste0("Tabula 10X/", "annotations_droplet.csv"), sep = ",", header = TRUE))
annotations <- annotations[annotations$cell %in% as.list(cellbarcodes), ]

tabulacol <- colnames(tabulamat)
tabulamat <- tabulamat[, tabulacol %in% annotations$cell]
tabulamat <- tabulamat[, order(colnames(tabulamat))]
annotations <- annotations[order(annotations$cell), ]
rm(tabulacol)

tabula10X <- CreateSeuratObject(counts = tabulamat, project = "tabula10X")
tabula10X@meta.data$tissue <- annotations$tissue
tabula10X@meta.data$mouse_sex <- annotations$mouse.sex
tabula10X@meta.data$celltype <- annotations$cell_ontology_class
tabula10X@meta.data$celltypeid <- annotations$cell_ontology_id
tabula10X@meta.data$noclusters <- 0

saveRDS(tabula10X, "Tabula 10X/tabula10X.rds")
