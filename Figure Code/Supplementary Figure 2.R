library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(scran)
library(scRNAseq)
library(scater)
library(BBmisc)
library(cowplot)
library(tidyverse)

# Panel A - B
mesc_ercc <- readRDS("Kolodziejczyk2015/mesc_ercc")
mesc_log <- readRDS("Kolodziejczyk2015/mesc_log")
mesc_log @meta.data$nCount_logRNA <- colSums(mesc_log @assays$RNA @data)
mesc_ercc @meta.data$celltype <- Idents(mesc_ercc)
mesc_log @meta.data$celltype <- Idents(mesc_log)

par(pty = "s")
dimcols <- c("PRGn")

DotPlot(mesc_ercc, features = c(
  "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Top2b",
  "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
  "Rpl9", "Gapdh", "Actb", "Sox2", "Pou5f1"
), cols = "PRGn") +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.35, axis.text = element_text(size = 13),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

DotPlot(mesc_log, features = c(
  "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Top2b",
  "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
  "Rpl9", "Gapdh", "Actb", "Sox2", "Pou5f1"
), cols = "PRGn") +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.35, axis.text = element_text(size = 13),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 10, unit = "pt")
  )

# Panel C - D
gonads_log <- readRDS("Niu2020/gonads_glo.rds")
gonads <- readRDS("Niu2020/gonads_abs.rds")

# Dotplot
Idents(gonads) <- "timepoint"
DotPlot(gonads,
  features = c(
    "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Mycn", "Mycl", "Max", "Top2b",
    "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
    "Rpl9", "Gapdh", "Actb", "Pou5f1"
  ),
  split.by = "germ", cols = c("PRGn")
) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.5, axis.text = element_text(size = 13),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

Idents(gonads_log) <- "timepoint"
DotPlot(gonads_log,
  features = c(
    "Chd1", "Hira", "Ino80", "Kat5", "Myc", "Mycn", "Mycl", "Max", "Top2b",
    "Atm", "Parp1", "Rps2", "Rpl3", "Rpl4", "Rpl6",
    "Rpl9", "Gapdh", "Actb", "Pou5f1"
  ),
  split.by = "germ", cols = c("PRGn")
) +
  ylab("") + xlab("") + RotatedAxis() +
  theme(
    aspect.ratio = 0.5, axis.text = element_text(size = 13),
    plot.title = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )
