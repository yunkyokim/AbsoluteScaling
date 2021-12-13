library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)
library(viridis)
library(stats)
library(ggridges)

#setwd("current_dir/")
tabulafacs = readRDS("tabulafacs.rds")
tabula10X = readRDS("tabula10X.rds")

#Tabula 10X
figcelltype = as.data.frame(tabula10X@meta.data$celltype)
colnames(figcelltype) = "CellType"
figcelltype[, "Transcripts"] = tabula10X@meta.data$nCount_RNA
figcelltype[, "Organ"] = tabula10X@meta.data$tissue

#Color Coordination
p = VlnPlot(
  tabula10X,
  features = c("nCount_RNA"),
  group.by = 'tissue',
  pt.size = 0,
  sort = FALSE
) + NoLegend() + theme(plot.margin = unit(c(0, 0, 0, 4), "cm"))
pbuild = ggplot2::ggplot_build(p)
pdata = pbuild$data[[1]]
cols = as.data.frame(unique(tabula10X@meta.data$tissue))
colnames(cols) = "Organs"
cols$colors = cols$Organs
cols = cols[order(cols$Organs), ]
cols[, "colors"] = as.data.frame(unique(pdata$fill))

figcelltype$colors = cols$colors[match(figcelltype$Organ, cols$Organs)]
figcelltype2 = figcelltype[(figcelltype$CellType != ""), ]
figcelltype2[, "organcell"] = paste(figcelltype2$CellType, "/", figcelltype2$Organ, sep = "")

filtercell <-
  as.data.frame(paste(figcelltype2$CellType, "/", figcelltype2$Organ, sep = ""))
colnames(filtercell) = "organcell"
filtercell$Transcripts = figcelltype2$Transcripts
filtercell = aggregate(. ~ organcell, filtercell, mean)
filtercell[, "CellType"] = sapply(X = strsplit(filtercell$organcell, split = "/"),
                                  FUN = "[",
                                  1)
filtercell = filtercell[order(filtercell$Transcripts), ]
filtercell[, "Duplicated"] = duplicated(filtercell$CellType)

filtercell = filtercell[filtercell$Duplicated == TRUE, ]
figcelltype2 = figcelltype2[!(figcelltype2$organcell %in% filtercell$organcell), ]
figcelltype$Organ = gsub("_", " ", figcelltype$Organ)

#Panel A, 450 x 450
ggplot(figcelltype, aes(
  x = reorder(Organ,-Transcripts, median),
  y = Transcripts,
  fill = Organ
)) +  geom_violin(trim = FALSE, scale = "width") + theme_classic() +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 5,
    colour = "grey27",
    shape = 95
  ) +
  NoLegend() + theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12)
  ) + labs(title = "", y = "Transcripts", x = "")

#Panel B, 650 x 1100
ggplot(figcelltype2, aes(
  y = reorder(CellType, Transcripts),
  x = Transcripts,
  fill = Organ
)) +  geom_density_ridges(scale = 2) + scale_y_discrete(position = "right") +
  theme_classic() + theme() + NoLegend() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(color = "black")) + labs(title = "", y = "", x = "Transcripts")

#Stats
organstats = data.frame(
  transcripts = tabula10X@meta.data$nCount_RNA,
  Identity = tabula10X@meta.data$tissue
)
organstats = setDT(organstats)[, list(
  Count = length(transcripts),
  Min = min(transcripts),
  Max = max(transcripts),
  Median = median(transcripts),
  IQR = IQR(transcripts)
), by = list(Identity)]

celltypestats = data.frame(
  transcripts = tabula10X@meta.data$nCount_RNA,
  Identity = tabula10X@meta.data$celltype
)
celltypestats = setDT(celltypestats)[, list(
  Count = length(transcripts),
  Min = min(transcripts),
  Max = max(transcripts),
  Median = median(transcripts),
  IQR = IQR(transcripts)
), by = list(Identity)]
celltypestats = celltypestats[celltypestats$Count > 1, ]

write.csv(organstats, "tabula10X_organstats.csv")
write.csv(celltypestats, "tabula10X_celltypestats.csv")


#Global UMAPs
tabula_unnorm = readRDS("tabulaunnorm.rds")
tabula_unnorm@meta.data$logRNA = tabulafacs@meta.data$logRNA
tabula_unnorm@meta.data$facsRNA = tabulafacs@meta.data$nCount_RNA

tabula_unnorm = NormalizeData(tabula_unnorm, normalization.method  = "RC")
tabula_unnorm@assays$RNA@data = as.matrix(log2(tabula_unnorm@assays$RNA@counts +
                                                 1))
tabula_unnorm <-
  FindVariableFeatures(tabula_unnorm,
                       selection.method = "vst",
                       nfeatures = 2000)
all.genes <- rownames(tabula_unnorm)
tabula_unnorm <- ScaleData(tabula_unnorm, features = all.genes)
tabula_unnorm <-
  RunPCA(tabula_unnorm, features = VariableFeatures(object = tabula_unnorm))
tabula_unnorm <- FindNeighbors(tabula_unnorm, dims = 1:15)
tabula_unnorm <- FindClusters(tabula_unnorm, resolution = 0.4)
tabula_unnorm <-
  RunUMAP(tabula_unnorm, dims = 1:15, min.dist = 0.75)

#Panels C
FeaturePlot(tabula_unnorm,
            features = c('logRNA'),
            pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Tabula Muris (Smart-seq2)")

DimPlot(
  tabula_unnorm,
  reduction = "umap",
  group.by = "tissue",
  pt.size = 0.1
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Tabula Muris (Smart-seq2)")


#Panels D
tabula10X@meta.data$logRNA = log2(tabula10X@meta.data$nCount_RNA + 1)

tabula10X = NormalizeData(tabula10X, normalization.method  = "RC")
tabula10X@assays$RNA@data = as.matrix(log2(tabula10X@assays$RNA@counts +
                                             1))
tabula10X <-
  FindVariableFeatures(tabula10X, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tabula10X)
tabula10X <- ScaleData(tabula10X, features = all.genes)
tabula10X <-
  RunPCA(tabula10X, features = VariableFeatures(object = tabula10X))
tabula10X <- FindNeighbors(tabula10X, dims = 1:15)
tabula10X <- FindClusters(tabula10X, resolution = 0.4)
tabula10X <- RunUMAP(tabula10X, dims = 1:15, min.dist = 0.75)

FeaturePlot(tabula10X, features = c('logRNA'), pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Tabula Muris (10X)") #450 x 350

DimPlot(
  tabula10X,
  reduction = "umap",
  group.by = "tissue",
  pt.size = 0.1
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  labs(title = "Tabula Muris (10X)")
