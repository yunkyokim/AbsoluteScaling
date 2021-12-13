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

#setwd("current_dir/")
tabula_intestine = readRDS("intestine.rds")
spliced_ercc = readRDS("spliced_ercc.rds")
unspliced_ercc = readRDS("unspliced_ercc.rds")

tabula_intestine@meta.data$cellnames = colnames(tabula_intestine)
tabula_intestine = subset(tabula_intestine, subset = cellnames %in% colnames(spliced_ercc))

tabula_intestine[["spliced"]] <-
  CreateAssayObject(counts = spliced_ercc)
tabula_intestine[["unspliced"]] <-
  CreateAssayObject(counts = unspliced_ercc)
rm(list = ls()[!ls() %in% c("tabula_intestine")])

tabula_intestine@assays$RNA@data = as.matrix(log2(tabula_intestine@assays$RNA@data +
                                                    1))
tabula_intestine <-
  FindVariableFeatures(tabula_intestine,
                       selection.method = "vst",
                       nfeatures = 2000)
all.genes <- rownames(tabula_intestine)
tabula_intestine <-
  ScaleData(tabula_intestine, features = all.genes)
tabula_intestine <-
  RunPCA(tabula_intestine, features = VariableFeatures(object = tabula_intestine))
tabula_intestine <- FindNeighbors(tabula_intestine, dims = 1:15)
tabula_intestine <- FindClusters(tabula_intestine, resolution = 0.4)
tabula_intestine <-
  RunUMAP(tabula_intestine, dims = 1:15, min.dist = 0.75)

#Panels G, 500 x 300
FeaturePlot(tabula_intestine,
            features = c('nCount_RNA'),
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
  labs(title = "Spliced Transcripts")

FeaturePlot(tabula_intestine,
            features = c('nCount_unspliced'),
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
  labs(title = "Unspliced Transcripts")

FeaturePlot(tabula_intestine,
            features = c('Chd1'),
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
  labs(title = "Chd1")

FeaturePlot(tabula_intestine,
            features = c('Myc'),
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
  labs(title = "Myc")

#Panels H, 500 x 300
vircols <- viridis(10, direction = 1, option = "D")
FeaturePlot(
  tabula_intestine,
  features = c('GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I'),
  pt.size = 0.1,
  cols = vircols
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
  labs(title = "GO-TRANSCRIPTION-BY-RNA-POLYMERASE-I")

FeaturePlot(
  tabula_intestine,
  features = c('Hypertranscription'),
  pt.size = 0.1,
  cols = vircols
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
  labs(title = "Serum Hypertranscription")

FeaturePlot(
  tabula_intestine,
  features = c('GO-RIBOSOME-BIOGENESIS'),
  pt.size = 0.1,
  cols = vircols
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
  labs(title = "GO-RIBOSOME-BIOGENESIS")

FeaturePlot(
  tabula_intestine,
  features = c('GO-DNA-REPAIR'),
  pt.size = 0.1,
  cols = vircols
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
  labs(title = "GO-DNA-REPAIR")

#Panel A, 950 x 300
DimPlot(
  tabula_intestine ,
  reduction = "umap",
  pt.size = 0.1,
  group.by = "celltype"
) +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

tabulacols = read.csv("tabulacols.csv")
organcols = rep("#00BF7D", 50)

#Panel D, 650 x 1100
library(ggridges)
figcelltype = as.data.frame(tabula_intestine@meta.data$celltype)
colnames(figcelltype) = "CellType"
figcelltype[, "Transcripts"] = tabula_intestine@meta.data$logRNA

ggplot(figcelltype, aes(y = reorder(CellType, Transcripts, median), x =
                          Transcripts)) +  geom_density_ridges(scale = 2, fill = "#00BF7D") + scale_y_discrete(position = "right") +
  theme_classic() + theme() + NoLegend() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(color = "black")) + labs(title = "", y = "", x = "log2 Transcripts") #500 x 250

#Panel F, 880 x 700
Idents(tabula_intestine) = "celltype"
tabula_intestine = ReorderIdent(object = tabula_intestine, var = "nCount_RNA", afxn = median)
DotPlot(
  tabula_intestine,
  features = c(
    "Chd1",
    "Hira",
    "Ino80",
    "Kat5",
    "Myc",
    "Top2b",
    "Atm",
    "Parp1",
    "Rps2",
    "Rpl3",
    "Rpl4",
    "Rpl6",
    "Rpl9",
    "Gapdh",
    "Actb",
    "Lgr5",
    "Krt20",
    "Atoh1",
    "Chga"
  )
) +
  labs(title = "", y = "", x = "") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(size = 10.5)
  ) +
  NoLegend()

#Gene Curves
intestine_curve = as.data.frame(t(as.data.frame(tabula_intestine@assays$RNA@data)))
intestine_curve[, "celltype"] = tabula_intestine@meta.data$celltype
intestine_curve = aggregate(. ~ celltype, intestine_curve, mean)
rownames(intestine_curve) = intestine_curve$celltype
intestine_curve$celltype = NULL
intestine_curve = as.data.frame(t(intestine_curve))
intestine_curve[, "total"] = rowSums(intestine_curve)
intestine_curve[, "genes"] = rownames(intestine_curve)
intestine_curve = intestine_curve[intestine_curve$total > 0, ]
intestine_curve = intestine_curve[order(intestine_curve$total), ]

library(tidyverse)
order.total = order(intestine_curve$total, decreasing = TRUE)
intestine_curve$rank <- NA
intestine_curve$rank[order.total] <- 1:nrow(intestine_curve)

intestine_tidy = intestine_curve %>%
  select(-c(total, genes)) %>%
  gather(key = "variable", value = "value",-rank)

#Panels E, 600 x 350
intestinecols = c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "#00BF7D")
ggplot(intestine_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 7,
      r = 0,
      b = 7,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Goblet Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = intestinecols) +
  NoLegend()

intestinecols = c("lightgrey", "lightgrey", "lightgrey", "#00BF7D", "lightgrey")
ggplot(intestine_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 7,
      r = 0,
      b = 7,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Epithelial Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = intestinecols) +
  NoLegend()

intestinecols = c("#00BF7D", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
ggplot(intestine_tidy, aes(x = rank, y = value)) +
  geom_smooth(aes(color = variable), size = 1.2, alpha = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 7,
      r = 0,
      b = 7,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "Enteroendocrine Cell", y = "Normalized Expression", x = "Ranked Genes") +
  scale_color_manual(values = intestinecols) +
  NoLegend()

#Lgr5 Annotation
Lgr5 = as.data.frame(tabula_intestine@assays$RNA@counts["Lgr5", ])
colnames(Lgr5) = "Lgr5"
for (row in 1:nrow(Lgr5)) {
  expr <- Lgr5[row, "Lgr5"]
  if (expr > 2.5) {
    Lgr5[row, "status"] = "Lgr5+"
  } else {
    Lgr5[row, "status"] = "Lgr5-"
  }
}
tabula_intestine@meta.data$lgr5 = Lgr5$status

#Panel B
DimPlot(
  tabula_intestine ,
  reduction = "umap",
  pt.size = 0.1,
  group.by = "lgr5",
  cols = c("lightgrey", "#00BF7D")
) +
  theme(aspect.ratio = 1, legend.text = element_text(size = 11)) +
  labs(title = "")

#Panel C
VlnPlot(
  tabula_intestine,
  features = c("nCount_RNA"),
  group.by = "lgr5" ,
  cols = c("lightgrey", "#00BF7D"),
  ncol = 1,
  pt.size = 0
) +
  ylab("Transcript Count") + xlab("") + NoLegend() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) # 600 x 350