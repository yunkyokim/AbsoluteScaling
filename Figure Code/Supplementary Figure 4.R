library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)
library(viridis)
library(stats)
library(ggridges)

tabula10X <- readRDS("Tabula 10X/tabula10X.rds")

# Tabula 10X
figcelltype <- as.data.frame(tabula10X@meta.data$celltype)
colnames(figcelltype) <- "CellType"
figcelltype[, "Transcripts"] <- tabula10X@meta.data$nCount_RNA
figcelltype[, "Organ"] <- tabula10X@meta.data$tissue

# Color Coordination
p <- VlnPlot(
  tabula10X,
  features = c("nCount_RNA"),
  group.by = "tissue",
  pt.size = 0,
  sort = FALSE
) + NoLegend() + theme(plot.margin = unit(c(0, 0, 0, 4), "cm"))
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]
cols <- as.data.frame(unique(tabula10X@meta.data$tissue))
colnames(cols) <- "Organs"
cols$colors <- cols$Organs
cols <- cols[order(cols$Organs), ]
cols[, "colors"] <- as.data.frame(unique(pdata$fill))
write.csv(cols, "Tabula 10X/tabula10Xcols.csv")

figcelltype$colors <- cols$colors[match(figcelltype$Organ, cols$Organs)]
figcelltype2 <- figcelltype[(figcelltype$CellType != ""), ]
figcelltype2[, "organcell"] <- paste(figcelltype2$CellType, "/", figcelltype2$Organ, sep = "")

filtercell <-
  as.data.frame(paste(figcelltype2$CellType, "/", figcelltype2$Organ, sep = ""))
colnames(filtercell) <- "organcell"
filtercell$Transcripts <- figcelltype2$Transcripts
filtercell <- aggregate(. ~ organcell, filtercell, mean)
filtercell[, "CellType"] <- sapply(
  X = strsplit(filtercell$organcell, split = "/"),
  FUN = "[",
  1
)
filtercell <- filtercell[order(filtercell$Transcripts), ]
filtercell[, "Duplicated"] <- duplicated(filtercell$CellType)

filtercell <- filtercell[filtercell$Duplicated == TRUE, ]
figcelltype2 <- figcelltype2[!(figcelltype2$organcell %in% filtercell$organcell), ]
figcelltype$Organ <- gsub("_", " ", figcelltype$Organ)

# Panel A,
ggplot(figcelltype, aes(
  x = reorder(Organ, -Transcripts, median),
  y = Transcripts,
  fill = Organ
)) +
  geom_violin(trim = FALSE, scale = "width") +
  theme_classic() +
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 5,
    colour = "grey27",
    shape = 95
  ) +
  NoLegend() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "", y = "Transcripts", x = "")

# Panel B,
ggplot(figcelltype2, aes(
  y = reorder(CellType, Transcripts),
  x = Transcripts,
  fill = Organ
)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme() +
  NoLegend() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_text(color = "black")
  ) +
  labs(title = "", y = "", x = "Transcripts")

# Stats
organstats <- data.frame(
  transcripts = tabula10X@meta.data$nCount_RNA,
  Identity = tabula10X@meta.data$tissue
)
organstats <- setDT(organstats)[, list(
  Count = length(transcripts),
  Min = min(transcripts),
  Max = max(transcripts),
  Median = median(transcripts),
  IQR = IQR(transcripts)
), by = list(Identity)]

celltypestats <- data.frame(
  transcripts = tabula10X@meta.data$nCount_RNA,
  Identity = tabula10X@meta.data$celltype
)
celltypestats <- setDT(celltypestats)[, list(
  Count = length(transcripts),
  Min = min(transcripts),
  Max = max(transcripts),
  Median = median(transcripts),
  IQR = IQR(transcripts)
), by = list(Identity)]
celltypestats <- celltypestats[celltypestats$Count > 1, ]

write.csv(organstats, "Tabula 10X/tabula10X_organstats.csv")
write.csv(celltypestats, "Tabula 10X/tabula10X_celltypestats.csv")

#Panel C
tabulafacs <- readRDS("Tabula FACS/tabulafacsabs.rds")
tabulafacs@meta.data$logRNA <- log2(tabulafacs@meta.data$nCount_RNA + 1)
tabula_unnorm <- readRDS("Tabula FACS/tabulaunnorm.rds") # No normalization data
tabula_unnorm@meta.data$logRNA <- tabulafacs@meta.data$logRNA
tabula_unnorm@meta.data$facsRNA <- tabulafacs@meta.data$nCount_RNA

autoprocess <- function(organ, tabula) {
  organ.seurat <- subset(tabula, subset = tissue == organ)
  organ.seurat <- NormalizeData(organ.seurat, normalization.method = "RC")
  organ.seurat@assays$RNA@data <- as.matrix(log2(organ.seurat@assays$RNA@data +
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
organplot <- function(dataset, title, min, max) {
  FeaturePlot(dataset, features = c("logRNA"), pt.size = 0.1) +
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
    scale_color_gradientn(
      colours = c("lightgrey", "blue"),
      limits = c(min, max)
    ) +
    NoLegend()
}

#FACS Plots
organlist <- sort(unique(tabula_unnorm@meta.data$tissue))
plotnumber <- 1
for (organ in organlist) {
  title <- gsub("_", " ", organ)
  processed.seurat <- autoprocess(organ, tabula_unnorm)
  eval(parse(
    text = paste(
      "p",
      plotnumber,
      "= organplot(processed.seurat, title, 17, 23)",
      sep = ""
    )
  ))
  plotnumber <- plotnumber + 1
}
library(egg)
gridpanels <- paste0("p", as.list(1:(plotnumber - 1)), collapse = ",")
eval(parse(text = paste(
  "grid.arrange(", gridpanels, ",ncol = 5)",
  sep = ""
)))



tabula10X <- readRDS("Tabula 10X/tabula10X.rds")
tabula10X@meta.data$logRNA <- log2(tabula10X@meta.data$nCount_RNA + 1)

#10X Plots
organlist <- sort(unique(tabula10X@meta.data$tissue))
plotnumber <- 1
for (organ in organlist) {
  title <- gsub("_", " ", organ)
  processed.seurat <- autoprocess(organ, tabula10X)
  eval(parse(
    text = paste(
      "p",
      plotnumber,
      "= organplot(processed.seurat, title, 10, 16)",
      sep = ""
    )
  ))
  plotnumber <- plotnumber + 1
}
library(egg)
gridpanels <- paste0("p", as.list(1:(plotnumber - 1)), collapse = ",")
eval(parse(text = paste(
  "grid.arrange(", gridpanels, ",ncol = 5)",
  sep = ""
)))


# Plots for legend/scales
FeaturePlot(processed.seurat,
            features = c("logRNA"),
            pt.size = 0.1
) +
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
    colours = c("lightgrey", "blue"),
    limits = c(10, 16),
    breaks = waiver(),
    n.breaks = 3
  )

FeaturePlot(processed.seurat,
            features = c("logRNA"),
            pt.size = 0.1
) +
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
    colours = c("lightgrey", "blue"),
    limits = c(17, 23),
    n.breaks = 4
  )
