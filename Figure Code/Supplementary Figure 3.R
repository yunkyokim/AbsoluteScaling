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

tabulafacs <- readRDS("Tabula FACS/tabulafacsabs.rds")
tabula_unnorm <- readRDS("Tabula FACS/tabulaunnorm.rds")

tabulafacs@meta.data$noclusters <- 0
VlnPlot(
  tabulafacs,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "noclusters",
  ncol = 3,
  pt.size = 0.1
)

# Figure Panels
# Panel E, 600 x 350
VlnPlot(
  tabulafacs,
  features = c("nCount_RNA"),
  cols = "maroon",
  ncol = 1,
  pt.size = 0,
  group.by = "noclusters"
) +
  labs(y = "Transcript Count", x = "Post Spike-in Normalization", title = "Tabula FACS") + NoLegend() +
  theme(
    aspect.ratio = 1.2,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_blank()
  )

# Panel C
VlnPlot(
  tabula_unnorm,
  features = c("nCount_RNA"),
  cols = "maroon",
  ncol = 1,
  pt.size = 0,
  group.by = "noclusters"
) +
  labs(y = "Transcript Count", x = "Pre Spike-in Normalization", title = "Tabula FACS") + NoLegend() +
  theme(
    aspect.ratio = 1.2,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_blank()
  )

# Panel D
VlnPlot(
  tabulafacs,
  features = c("nFeature_RNA"),
  cols = "maroon",
  ncol = 1,
  pt.size = 0,
  group.by = "noclusters"
) +
  labs(y = "Feature Count", x = "", title = "Tabula FACS") + NoLegend() +
  theme(
    aspect.ratio = 1.2,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_blank()
  )


# ERCC Correlation Matrix
genelist <- rownames(tabula_unnorm)
ercclist <- (genelist[grep(c("ERCC-"), genelist)])
erccmat <- as.data.frame(tabula_unnorm@assays$RNA@counts)
erccmat <- erccmat[ercclist, ]
erccmat <- as.data.frame(t(erccmat))
erccmat$tissue <- tabula_unnorm@meta.data$tissue
erccmat <- aggregate(. ~ tissue, erccmat, sum)
ercc_cor <- erccmat
rownames(ercc_cor) <- ercc_cor$tissue
ercc_cor$tissue <- NULL
ercc_cor <- as.data.frame(t(ercc_cor))
colnames(ercc_cor) <- gsub("_", " ", colnames(ercc_cor))
cormat <- round(cor(ercc_cor), 6)

library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

ggplot(data = melted_cormat, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile()
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri


melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Panel G
library(RColorBrewer)
newcol <- colorRampPalette((brewer.pal(9, "RdPu")))
puor <- newcol(200)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = puor[1:170], name = "Pearson\nCorrelation") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 11,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.3, 0.5),
    legend.direction = "vertical",
    legend.text = element_text(size = 10),
    aspect.ratio = 1
  ) +
  labs(title = "ERCC Spike-in Expression") +
  guides(
    fill = guide_colorbar(
      barwidth = 2,
      barheight = 7,
      title.position = "top",
      title.hjust = 0.5,
      title.vjust = 3
    )
  )


# Panel F
VlnPlot(
  tabula_unnorm,
  features = c("ERCC"),
  cols = "maroon",
  ncol = 1,
  pt.size = 0,
  group.by = "noclusters"
) +
  labs(y = "ERCC Percentage", x = "", title = "Tabula FACS") + NoLegend() +
  theme(
    aspect.ratio = 1.2,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_blank()
  )

ercc_sum <- as.data.frame(tabula_unnorm@assays$RNA@counts[ercclist, ])
ercc_sum <- colSums(ercc_sum)
tabula_unnorm@meta.data$ERCC_SUM <- ercc_sum

library(MASS)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

ercc_scatter <- as.data.frame(tabula_unnorm@meta.data$nCount_RNA)
colnames(ercc_scatter) <- "RNA_Counts"
ercc_scatter[, "ERCC"] <- tabula_unnorm@meta.data$ERCC_SUM
ercc_scatter$ERCCnorm <- BBmisc::normalize(
  ercc_scatter$ERCC,
  method = "range",
  range = c(0, 10),
  margin = 2L
)
ercc_scatter$RNA_Countsnorm <- BBmisc::normalize(
  ercc_scatter$RNA_Counts,
  method = "range",
  range = c(0, 10),
  margin = 2L
)

ercc_scatter$density <-
  get_density(ercc_scatter$ERCC, ercc_scatter$RNA_Counts, n = 300)

# Panel H
ggplot(ercc_scatter, aes(x = ERCC, y = RNA_Counts)) +
  theme_classic() +
  geom_point(aes(ERCC, RNA_Counts, color = density),
    size = 0.01,
    alpha = 1
  ) +
  theme(
    aspect.ratio = 1,
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
  labs(title = "", y = "Total Transcripts", x = "ERCC Transcripts") +
  lims(x = c(300, 350000), y = c(50000, 6000000)) +
  scale_color_viridis(name = "Density")
NoLegend()

ercccontrols <- read.delim2("cms_095046.txt") # ERCC concentration data from Thermo Fisher
ercccontrols <- ercccontrols[order(ercccontrols$ERCC.ID), ]
ercccontrols$conc <- ercccontrols$concentration.in.Mix.1..attomoles.ul.

ercc_conc <- erccmat
rownames(ercc_conc) <- erccmat$tissue
ercc_conc$tissue <- NULL
ercc_conc <- as.data.frame(t(ercc_conc))
ercc_conc$ERCC <- rowSums(ercc_conc)
ercc_conc$ERCCnames <- rownames(ercc_conc)
ercc_conc$conc <- as.numeric(ercccontrols$conc)

# Panel I
ggplot(ercc_conc, aes(x = conc, y = ERCC)) +
  theme_classic() +
  geom_point() +
  geom_smooth(
    method = lm,
    colour = "violetred3",
    size = 1
  ) +
  theme(
    aspect.ratio = 1,
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
  labs(title = "", y = "Total Log2 ERCC Transcripts", x = "Expected Log2 ERCC Transcripts (attomoles/uL)") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
  scale_x_continuous(trans = "log2", labels = scales::scientific) +
  scale_y_continuous(trans = "log2", labels = scales::scientific) +
  NoLegend()

range(ercc_conc$conc)
range(ercc_conc$ERCC)


# Tabula 10X Comparison
tabula10X <- readRDS("Tabula 10X/tabula10X.rds")

tabula10X_RNA <- as.data.frame(tabula10X@meta.data$nCount_RNA, row.names = tabula10X@meta.data$tissue)
tabula10X_RNA[, "tissue"] <- rownames(tabula10X_RNA)
tabula10X_RNA <- aggregate(. ~ tissue, tabula10X_RNA, median)
names(tabula10X_RNA)[names(tabula10X_RNA) == "tabula10X@meta.data$nCount_RNA"] <-
  "UMI_Counts"
rownames(tabula10X_RNA) <- tabula10X_RNA[, "tissue"]

tabulafacs_RNA <- as.data.frame(tabulafacs@meta.data$nCount_RNA,
  row.names = tabulafacs@meta.data$tissue
)
tabulafacs_RNA[, "tissue"] <- rownames(tabulafacs_RNA)
tabulafacs_RNA <- aggregate(. ~ tissue, tabulafacs_RNA, median)
names(tabulafacs_RNA)[names(tabulafacs_RNA) == "tabulafacs@meta.data$nCount_RNA"] <-
  "ERCC_Counts"
rownames(tabulafacs_RNA) <- tabulafacs_RNA[, "tissue"]
tabulafacs_RNA$tissue <- gsub("_", " ", tabulafacs_RNA$tissue)

tabulafacs_RNA <- subset(tabulafacs_RNA, tissue %in% tabula10X_RNA$tissue)
tabula10X_RNA <- subset(tabula10X_RNA, tissue %in% tabulafacs_RNA$tissue)
tabulafacs_RNA[, "UMI_Counts"] <- tabula10X_RNA$UMI_Counts

plot <- ggplot(tabulafacs_RNA, aes(x = ERCC_Counts, y = UMI_Counts)) +
  theme_classic() +
  geom_point(aes(color = rownames(tabulafacs_RNA)), size = 3) +
  geom_smooth(
    color = "black",
    fill = "grey90",
    linetype = "dashed",
    method = "lm",
    se = FALSE
  ) +
  geom_text_repel(
    aes(
      label = rownames(tabulafacs_RNA),
      color = rownames(tabulafacs_RNA)
    ),
    size = 3,
    point.padding = 0.5,
    segment.alpha = 0
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
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
  labs(title = "", y = "Median Transcripts (Tabula 10X)", x = "Median Transcripts (Tabula FACS)") +
  stat_cor(method = "pearson", size = 4) +
  scale_x_continuous(trans = "log2", labels = scales::scientific) +
  scale_y_continuous(trans = "log2", labels = scales::scientific) +
  NoLegend()
plot

# Panel K
g <- ggplot_build(plot)
organcolors <- as.data.frame(g$data[[1]]$colour, row.names = g$data[[3]]$label)
organcolors <- as.data.frame(t(organcolors))

# Panel A
VlnPlot(
  tabula10X,
  features = c("nCount_RNA"),
  cols = "lightblue",
  ncol = 1,
  pt.size = 0,
  group.by = "noclusters"
) +
  labs(y = "Transcript Count", x = "", title = "Tabula 10X") + NoLegend() +
  theme(
    aspect.ratio = 1.2,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_blank()
  ) # 600 x 350

# Panel B
VlnPlot(
  tabula10X,
  features = c("nFeature_RNA"),
  cols = "lightblue",
  ncol = 1,
  pt.size = 0,
  group.by = "noclusters"
) +
  labs(y = "Feature Count", x = "", title = "Tabula 10X") + NoLegend() +
  theme(
    aspect.ratio = 1.2,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_blank()
  ) # 600 x 350


tabulafacs@assays$RNA@data <- as.matrix(log2(tabulafacs@assays$RNA@counts +
  1))
tabulafacs_GENE <- as.data.frame(tabulafacs@assays$RNA@data)
tabulafacs_GENE <- as.data.frame(t(tabulafacs_GENE))
tabulafacs_GENE[, "tissue"] <- tabulafacs@meta.data$tissue
tabulafacs_GENE <- aggregate(. ~ tissue, tabulafacs_GENE, median)
tabulafacs_GENE$tissue <- paste0("FACS_", tabulafacs_GENE$tissue)

tabula10X@assays$RNA@data <- as.matrix(log2(tabula10X@assays$RNA@counts +
  1))
tabula10X_GENE <- as.data.frame(tabula10X@assays$RNA@data)
tabula10X_GENE <- as.data.frame(t(tabula10X_GENE))
tabula10X_GENE[, "tissue"] <- tabula10X@meta.data$tissue
tabula10X_GENE <- aggregate(. ~ tissue, tabula10X_GENE, median)
tabula10X_GENE$tissue <- paste0("10X", tabula10X_GENE$tissue)

tabula10X_GENE$tissue <- gsub(" ", "_", tabula10X_GENE$tissue)
rownames(tabula10X_GENE) <- tabula10X_GENE$tissue

tabulafacs_GENE <- tabulafacs_GENE[, intersect(colnames(tabula10X_GENE), colnames(tabulafacs_GENE))]
tabula10X_GENE <- tabula10X_GENE[, intersect(colnames(tabula10X_GENE), colnames(tabulafacs_GENE))]
rownames(tabulafacs_GENE) <- tabulafacs_GENE$tissue
rownames(tabula10X_GENE) <- tabula10X_GENE$tissue
tabulafacs_GENE[rownames(tabula10X_GENE), ] <- tabula10X_GENE[rownames(tabula10X_GENE), ]
tabulafacs_GENE$tissue <- NULL
tabulafacs_GENE <- t(tabulafacs_GENE)
tabulafacs_GENE <- as.data.frame(tabulafacs_GENE)

gene_comparison_plot <- function(organ, title) {
  eval(parse(
    text = paste(
      "ggplot(subset(tabulafacs_GENE, `FACS_",
      organ,
      "` + `10X",
      organ,
      "`>0), aes(x = `FACS_",
      organ,
      "`, y = `10X",
      organ,
      "`)) +
    geom_point(color=organcolors$",
      organ,
      ',size = 1) +
    geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method="lm")+
    theme_classic() +
    theme(text=element_text(size=12),
    axis.text=element_text(size=10, colour = "black"),
    plot.title=element_text(size=11, hjust = 0.5),
    axis.title=element_text(size=12),
    aspect.ratio = 1) +
    labs(title = "',
      title,
      '", y = "Tabula 10X", x = "Tabula FACS")+
    stat_cor(method = "pearson", size = 3.5)+
    NoLegend()',
      sep = ""
    )
  ))
}

# Panel J
plotnumber <- 1
for (i in colnames(organcolors)) {
  title <- gsub("_", " ", i)
  eval(parse(
    text = paste("p", plotnumber, "= gene_comparison_plot(i,title)", sep = "")
  ))
  plotnumber <- plotnumber + 1
}
library(egg)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, nrow = 3)
