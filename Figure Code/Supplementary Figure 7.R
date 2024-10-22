library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(BBmisc)
library(viridis)

tabulafacs <- readRDS("Tabula FACS/tabulafacs.rds")
tabulafacs@meta.data$logRNA <- log2(tabulafacs@meta.data$nCount_RNA + 1)
tabula_unnorm <- readRDS("Tabula FACS/tabulaunnorm.rds")
tabula_unnorm@meta.data$logRNA <- tabulafacs@meta.data$logRNA
tabula_unnorm@meta.data$facsRNA <- tabulafacs@meta.data$nCount_RNA

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)

# Cell Cycle Scoring
tabula_unnorm <-
  NormalizeData(tabula_unnorm, normalization.method = "RC")
tabula_unnorm@assays$RNA@data <- as.matrix(log2(tabula_unnorm@assays$RNA@data +
  1))
tabula_unnorm <-
  CellCycleScoring(
    tabula_unnorm,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE
  )

phase_prop <- as.data.frame(tabula_unnorm@meta.data$celltype)
colnames(phase_prop) <- "celltype"
phase_prop$phase <- tabula_unnorm@meta.data$Phase
phase_prop$phase <- tabula_unnorm@meta.data$Phase

phase_prop2 <- as_tibble(phase_prop)

prop <- phase_prop2 %>% dplyr::count(celltype)
prop <- phase_prop2 %>%
  tidyr::gather(celltype, phase) %>%
  dplyr::group_by(celltype, phase) %>%
  dplyr::count() %>%
  tidyr::spread(celltype, n)
prop <- as.data.frame(prop)
rownames(prop) <- prop$phase
prop$phase <- NULL
prop <- as.data.frame(t(prop))

prop[is.na(prop)] <- 0
prop <- prop / rowSums(prop)
prop <- prop[-c(1), ]
prop$celltype <- rownames(prop)
prop$cycle <- prop$S + prop$G2M

prop$reordered <- reorder(prop$celltype, prop$cycle)
propmelt <-
  reshape2::melt(prop[, c("celltype", "G1", "G2M", "S", "reordered")], id.vars = c("celltype", "reordered"))
cols <- c(
  G1 = "gold2", G2M = "violetred3",
  S = "seagreen3"
)

# Panel A
ggplot(propmelt, aes(x = reordered, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 10),
    aspect.ratio = 2.5
  ) +
  scale_fill_manual(values = cols) +
  rotate() +
  labs(title = "", y = "Phase Percentage", x = "") +
  NoLegend()

cell_transcripts <- as.data.frame(tabula_unnorm@meta.data$facsRNA)
colnames(cell_transcripts) <- "transcripts"
cell_transcripts$celltype <- tabula_unnorm@meta.data$celltype
cell_transcripts <- aggregate(. ~ celltype, cell_transcripts, median)
cell_transcripts <- cell_transcripts[-1, ]
library(zoo)
prop[, "transcripts"] <- cell_transcripts$transcripts[match(cell_transcripts$celltype, prop$celltype)]
prop[, "transcript_roll"] <- zoo::rollmean(prop$transcripts, k = 7, fill = "extend")

# Panel B
ggplot(prop, aes(x = reordered, y = transcripts)) +
  geom_bar(stat = "identity", fill = "wheat3") +
  geom_line(aes(x = reordered, y = transcript_roll, group = 1),
    colour = "turquoise",
    size = 1
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 10),
    aspect.ratio = 2.5
  ) +
  rotate() +
  labs(title = "", y = "Transcripts", x = "")

prop <- data.frame(
  "s" = tabula_unnorm@meta.data$S.Score,
  "g2m" = tabula_unnorm@meta.data$G2M.Score,
  "phase" = tabula_unnorm@meta.data$Phase,
  "transcripts" = tabula_unnorm@meta.data$facsRNA
)
prop$cycle <- prop$s + prop$g2m
prop$celltype <- tabula_unnorm@meta.data$celltype
prop$id <- colnames(tabula_unnorm@assays$RNA@data)

library(MASS)
library(ggplot2)
library(viridis)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
prop$density <- get_density(prop$cycle, prop$transcripts, n = 200)

# Panel C
ggplot(prop, aes(x = cycle, y = transcripts, color = density)) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 10),
    aspect.ratio = 0.4
  ) +
  labs(
    title = "",
    y = "Log2 Transcripts",
    x = "Aggregate G2M/S Phase Score",
    color = "Density"
  )
