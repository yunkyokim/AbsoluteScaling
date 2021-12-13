library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)

#setwd("current_dir/")

#Data import
#Import normalized Seurat objects
gonads_log = readRDS("gonads_glo.rds")
gonads = readRDS("gonads_abs.rds")
Idents(gonads) = "timepoint"
Idents(gonads_log) = "timepoint"

#Panel C, 950 x 400
DotPlot(
  gonads_log,
  features = c(
    "Chd1",
    "Hira",
    "Ino80",
    "Kat5",
    "Myc",
    "Mycn",
    "Mycl",
    "Max",
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
    "Pou5f1"
  ),
  split.by = "germ",
  cols = c("blue", "blue")
) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.5,
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

#Panel D, 950 x 400
DotPlot(
  gonads,
  features = c(
    "Chd1",
    "Hira",
    "Ino80",
    "Kat5",
    "Myc",
    "Mycn",
    "Mycl",
    "Max",
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
    "Pou5f1"
  ),
  split.by = "germ",
  cols = c("blue", "blue")
) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.5,
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