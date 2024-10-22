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
tabula10X <- readRDS("Tabula 10X/tabula10X.rds")

# Tabula FACS
figcelltype <- as.data.frame(tabulafacs@meta.data$celltype)
colnames(figcelltype) <- "CellType"
figcelltype[, "Transcripts"] <- tabulafacs@meta.data$nCount_RNA
figcelltype[, "Organ"] <- tabulafacs@meta.data$tissue

# Color Coordination
p <- VlnPlot(
  tabulafacs,
  features = c("nCount_RNA"),
  group.by = "tissue",
  pt.size = 0,
  sort = FALSE
) + NoLegend() + theme(plot.margin = unit(c(0, 0, 0, 4), "cm"))
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]
cols <- as.data.frame(unique(tabulafacs@meta.data$tissue))
colnames(cols) <- "Organs"
cols$colors <- cols$Organs
cols <- cols[order(cols$Organs), ]
cols[, "colors"] <- as.data.frame(unique(pdata$fill))
write.csv(cols, "Tabula FACS/tabulafacscols.csv")

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

# Panel A
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
    axis.title = element_text(size = 12),
    aspect.ratio = 0.45
  ) +
  labs(title = "", y = "Transcripts", x = "")

# Panel B
library(ggridges)
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

# Organ Stats/Characterization
organstats <- data.frame(
  transcripts = tabulafacs@meta.data$nCount_RNA,
  Identity = tabulafacs@meta.data$tissue
)
organstats <- setDT(organstats)[, list(
  Count = length(transcripts),
  Min = min(transcripts),
  Max = max(transcripts),
  Median = median(transcripts),
  IQR = IQR(transcripts)
), by = list(Identity)]

celltypestats <- data.frame(
  transcripts = tabulafacs@meta.data$nCount_RNA,
  Identity = tabulafacs@meta.data$celltype
)
celltypestats <- setDT(celltypestats)[, list(
  Count = length(transcripts),
  Min = min(transcripts),
  Max = max(transcripts),
  Median = median(transcripts),
  IQR = IQR(transcripts)
), by = list(Identity)]
celltypestats <- celltypestats[celltypestats$Count > 1, ]

write.csv(organstats, "Tabula FACS/tabulafacs_organstats.csv")
write.csv(celltypestats, "Tabula Facs/tabulafacs_celltypestats.csv")

# Signature Analysis
signature.list <- ("c5.bp.v7.1.symbols.gmt") # from GSEA database
tabulafacsmat <- as.sparse(tabulafacs@assays$RNA@counts)
tabula.vision <- Vision(
  tabulafacsmat,
  signatures = signature.list,
  projection_methods = NULL,
  pool = FALSE,
  sig_gene_threshold = 0.05,
  sig_norm_method = "none"
)
tabula.vision <- calcSignatureScores(tabula.vision,
  sig_norm_method = "none",
  sig_gene_importance = FALSE
)

sigscores <- t(as.data.frame(getSignatureScores(tabula.vision)))
tabulafacs[["vision"]] <- CreateAssayObject(counts = sigscores)

serum_diff <- read.csv("Tabula FACS/mesc_serum_sig.csv", header = FALSE)
serum_diff[1, 1] <- "TagIn"
serum_list <- (serum_diff$V2)
names(serum_list) <- serum_diff$V1

hypertranscription <- createGeneSignature("Hypertranscription", sigData = serum_list)
tabula.vision <- Vision(
  tabulafacsmat,
  signatures = c(hypertranscription),
  projection_methods = NULL,
  pool = FALSE,
  sig_norm_method = "none"
)
tabula.vision <- calcSignatureScores(tabula.vision,
  sig_norm_method = "znorm_rows",
  sig_gene_importance = FALSE
)

sigscores <- t(as.data.frame(getSignatureScores(tabula.vision)))
tabulafacs[["hyper"]] <- CreateAssayObject(counts = sigscores)

SIGcomparison <- as.data.frame(t(tabulafacs@assays$vision@data))
sigscores <- as.data.frame(t(tabulafacs@assays$hyper@data))
SIGcomparison[, "Hypertranscription"] <- sigscores$Hypertranscription
SIGcomparison[, "RNA_Counts"] <- tabulafacs@meta.data$logRNA

saveRDS(tabulafacs, file = "Tabula FACS/tabulafacs_vision.rds")

sig <- colnames(SIGcomparison)
sig <- gsub(
  x = sig,
  pattern = "\\-",
  replacement = "_"
)
colnames(SIGcomparison) <- sig
rm(sig)

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

# Panel C
ggplot(
  SIGcomparison,
  aes(
    x = RNA_Counts,
    y = Hypertranscription,
    colour = get_density(RNA_Counts, Hypertranscription, n = 200)
  )
) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  geom_smooth(
    color = "black",
    fill = "grey90",
    linetype = "dashed",
    method = "lm",
    se = TRUE
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  NoLegend() +
  labs(y = "Serum Hypertranscription", x = "Log2 Transcripts")

# Panel D
ggplot(
  SIGcomparison,
  aes(
    x = RNA_Counts,
    y = GO_TRANSCRIPTION_BY_RNA_POLYMERASE_I,
    colour = get_density(RNA_Counts, GO_TRANSCRIPTION_BY_RNA_POLYMERASE_I, n = 200)
  )
) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  geom_smooth(
    color = "black",
    fill = "grey90",
    linetype = "dashed",
    method = "lm",
    se = TRUE
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  NoLegend() +
  labs(y = "GO_TRANSCRIPTION_BY_RNA_POLYMERASE_I", x = "Log2 Transcripts")

# Panel E
ggplot(
  SIGcomparison,
  aes(
    x = RNA_Counts,
    y = GO_RIBOSOME_BIOGENESIS,
    colour = get_density(RNA_Counts, GO_RIBOSOME_BIOGENESIS, n = 200)
  )
) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  geom_smooth(
    color = "black",
    fill = "grey90",
    linetype = "dashed",
    method = "lm",
    se = TRUE
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  NoLegend() +
  labs(y = "GO_RIBOSOME_BIOGENESIS", x = "Log2 Transcripts")

# Gene Correlation Comparison, Single-cell Level
tabulafacs@assays$RNA@data <- as.matrix(log2(tabulafacs@assays$RNA@counts +
  1))
GENEcomparison <- as.data.frame(t(tabulafacs@assays$RNA@data))
GENEcomparison[, "RNA_Counts"] <- tabulafacs@meta.data$logRNA

GENEtest <- data.frame()
for (sig in colnames(GENEcomparison)) {
  eval(parse(
    text = paste(
      "test = cor.test(GENEcomparison$RNA_Counts, GENEcomparison$`",
      sig,
      '`, method="spearman", exact = FALSE)',
      sep = ""
    )
  ))
  GENEtest[sig, "Statistic"] <- test$statistic
  GENEtest[sig, "P_Value"] <- test$p.value
  GENEtest[sig, "Estimate"] <- test$estimate
  GENEtest[sig, "Alternative"] <- test$alternative
  GENEtest[sig, "Method"] <- test$method
  print(sig)
}

GENEtest <- GENEtest[order(GENEtest$Statistic), ]
GENEtest$Correlation <- GENEtest$Estimate
GENEtest$Signature <- rownames(GENEtest)
GENEtest$logP_Value <- log(GENEtest$P_Value)
GENEtestPOS <- subset(GENEtest, P_Value <= 0.05)
write.table(GENEtestPOS, "Tabula FACS/geneSCtestpos.csv", sep = ",")

GENEtest$gene <- rownames(GENEtest)
genelist <- rownames(GENEtest)

ribosomepointlist <- c(genelist[grep(c("Rpl"), genelist)], genelist[grep(c("Rps"), genelist)])
GENEtest$gene <- rownames(GENEtest)
GENEtest <- GENEtest[order(-GENEtest$Statistic), ]

# Panel F
ggplot(GENEtest, aes(y = Estimate, x = reorder(gene, -Estimate))) +
  theme_classic() +
  geom_point(color = "grey", size = 2) +
  geom_label_repel(
    data = head(subset(GENEtest, gene %in% ribosomepointlist), n = 10),
    aes(label = gene),
    max.overlaps = Inf,
    box.padding = 0.5
  ) +
  geom_point(
    data = subset(GENEtest, gene %in% ribosomepointlist),
    color = "seagreen3",
    size = 2
  ) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    aspect.ratio = 1,
    axis.text.x = element_blank(),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(title = "", y = "Spearman Correlation", x = "Ranked Genes") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

# Gene Comparison, Single-Organ Level
# Panel G
GENEcomparison[, "Celltype"] <- tabulafacs@meta.data$celltype
GENEcomparison[, "Organ"] <- tabulafacs@meta.data$tissue
GENEcomparison[, "Organcell"] <- paste(GENEcomparison$Celltype, "/", GENEcomparison$Organ, sep = "")
GENEcomparison$Celltype <- NULL
GENEcomparison$Organ <- NULL

GENEcomparison <- aggregate(. ~ Organcell, GENEcomparison, mean)
GENEcomparison <- normalize(GENEcomparison, method = "range", range = c(0, 10), margin = 2L)
rownames(GENEcomparison) <- GENEcomparison$Organcell
GENEcomparison[, "Organcell"] <- NULL
GENEcomparison <- GENEcomparison[-c(1:3), ]

saveRDS(GENEcomparison, file = "Tabula FACS/tabulafacs_genecomparison.rds")

GENEcomparison[, "organcell"] <- rownames(GENEcomparison)
GENEcomparison[, "Organ"] <- sapply(X = strsplit(GENEcomparison$organcell, split = "/"), FUN = "[", 2)
GENEcomparison[, "Celltype"] <- sapply(X = strsplit(GENEcomparison$organcell, split = "/"), FUN = "[", 1)

GENEcomparison[, "colors"] <- cols$colors[match(GENEcomparison$Organ, cols$Organs)]
GENEcomparison$organcell <- paste(GENEcomparison$Celltype, " (", GENEcomparison$Organ, ")", sep = "")
rownames(GENEcomparison) <- GENEcomparison$organcell

ggplot(GENEcomparison, aes(x = RNA_Counts, y = Chd1)) +
  geom_point(aes(color = colors), size = 3) +
  geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm", se = TRUE) +
  geom_text_repel(aes(label = ifelse(RNA_Counts > 6 & Chd1 > 5.5, as.character(Celltype), ""), color = colors),
    size = 4, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_colour_identity() +
  labs(y = "Chd1", x = "Log2 Transcripts") + # 600 x 500
  NoLegend()

# Panel H
enrichrGO <- read.csv("Tabula FACS/enrichrGOBP.csv", fill = TRUE) # performed externally using highly correlated genes from GENEtest
enrichrGO$logp <- -log(enrichrGO$Adjusted.P.value)
enrichrGO$Name <- sapply(X = strsplit(enrichrGO$Term, split = " \\("), FUN = "[", 1)

ggplot(head(enrichrGO, 20), aes(x = logp, y = reorder(Name, logp)), ) +
  geom_bar(stat = "identity", fill = "turquoise4") +
  theme_classic() +
  NoLegend() +
  theme(
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 0,
      unit = "pt"
    ),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    aspect.ratio = 1.4,
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  xlim(0, 230) +
  labs(title = "GO BP Terms", y = "", x = "-log p-Value") +
  scale_y_discrete(
    label = function(x) {
      stringr::str_trunc(x, 45)
    }
  )
