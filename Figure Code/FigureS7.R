library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(scran)
library(BBmisc)

#setwd("current_dir/")
tabulafacs = readRDS("tabulafacs.rds")
tabulafacs@meta.data$logRNA = log2(tabulafacs@meta.data$nCount_RNA + 1)
tabula_unnorm = readRDS("tabulaunnorm.rds") #No normalization data
tabula_unnorm@meta.data$logRNA = tabulafacs@meta.data$logRNA
tabula_unnorm@meta.data$facsRNA = tabulafacs@meta.data$nCount_RNA

autoprocess = function(organ, tabula) {
  organ.seurat = subset(tabula, subset = tissue == organ)
  organ.seurat = NormalizeData(organ.seurat, normalization.method  = "RC")
  organ.seurat@assays$RNA@data = as.matrix(log2(organ.seurat@assays$RNA@counts +
                                                  1))
  organ.seurat <-
    FindVariableFeatures(organ.seurat,
                         selection.method = "vst",
                         nfeatures = 2000)
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
organplot = function(dataset, title, min, max) {
  FeaturePlot(dataset, features = c('logRNA'), pt.size = 0.1) +
    theme(
      aspect.ratio = 1,
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 10),
      plot.margin = margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1,
        unit = "pt"
      )
    ) +
    labs(title = title) +
    scale_color_gradientn(colours = c('lightgrey', 'blue'),
                          limits = c(min, max)) +
    NoLegend()
}

#Panel A, 1100 x 800
organlist = sort(unique(tabula_unnorm@meta.data$tissue))
plotnumber = 1
for (organ in organlist) {
  title = gsub('_', ' ', organ)
  processed.seurat = autoprocess(organ, tabula_unnorm)
  eval(parse(
    text = paste(
      "p",
      plotnumber,
      "= organplot(processed.seurat, title, 17, 23)",
      sep = ""
    )
  ))
  plotnumber = plotnumber + 1
}
library(egg)
gridpanels = paste0("p", as.list(1:(plotnumber - 1)), collapse = ",")
eval(parse(text = paste(
  "grid.arrange(", gridpanels, ",ncol = 5)", sep = ""
)))



tabula10X = readRDS("tabula10X.rds")
tabula10X@meta.data$logRNA = log2(tabula10X@meta.data$nCount_RNA + 1)

#Panel B, 1100 x 600
organlist = sort(unique(tabula10X@meta.data$tissue))
plotnumber = 1
for (organ in organlist) {
  title = gsub('_', ' ', organ)
  processed.seurat = autoprocess(organ, tabula10X)
  eval(parse(
    text = paste(
      "p",
      plotnumber,
      "= organplot(processed.seurat, title, 10, 16)",
      sep = ""
    )
  ))
  plotnumber = plotnumber + 1
}
library(egg)
gridpanels = paste0("p", as.list(1:(plotnumber - 1)), collapse = ",")
eval(parse(text = paste(
  "grid.arrange(", gridpanels, ",ncol = 5)", sep = ""
)))


#Plots for legend/scales
FeaturePlot(processed.seurat,
            features = c('logRNA'),
            pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  scale_color_gradientn(
    colours = c('lightgrey', 'blue'),
    limits = c(10, 16),
    breaks = waiver(),
    n.breaks = 3
  )

FeaturePlot(processed.seurat,
            features = c('logRNA'),
            pt.size = 0.1) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  scale_color_gradientn(
    colours = c('lightgrey', 'blue'),
    limits = c(17, 23),
    n.breaks = 4
  )
