library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)

#setwd("current_dir/")

#Data import
#Import normalized Seurat objects
mesc_ercc = readRDS("mesc_ercc")
mesc_log = readRDS("mesc_log")
mesc_log@meta.data$nCount_logRNA = colSums(mesc_log@assays$RNA@data)
mesc_ercc@meta.data$celltype = Idents(mesc_ercc)
mesc_log@meta.data$celltype = Idents(mesc_log)

#Figure Panels
par(pty = "s")
dimcols = c("mediumaquamarine", "gold3")

#Panel A, 800 x 300
DotPlot(
  mesc_log,
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
    "Sox2",
    "Pou5f1"
  )
) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.35,
    axis.text = element_text(size = 13),
    plot.title = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 10,
      unit = "pt"
    )
  )

#Panel B, 800 x 300
DotPlot(
  mesc_ercc,
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
    "Sox2",
    "Pou5f1"
  )
) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.35,
    axis.text = element_text(size = 13),
    plot.title = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  )
