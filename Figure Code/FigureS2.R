library(Seurat)
library(ggplot2)
library(BBmisc)

#setwd("current_dir/")

#Cell Hashings PBMC Analysis
hashing = readRDS("hashing.subset")

#Panel A Middle, 400 x 300
vlncol = c("hotpink4", "hotpink3")
VlnPlot(
  hashing,
  features = c("nCount_RNA"),
  group.by = 'orig.ident',
  cols = vlncol,
  ncol = 1,
  pt.size = 0
) + ylab("Transcripts") + xlab("") + NoLegend() + ggtitle("") +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 15,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#Fold-change calulation
hashing_transcripts = data.frame(
  transcripts = hashing@meta.data$nCount_RNA,
  identity = hashing@meta.data$orig.ident
)
hashing_transcripts = aggregate(. ~ identity, hashing_transcripts, median)
rownames(hashing_transcripts) = hashing_transcripts$identity
hashing_transcripts["Doublet", "transcripts"] / hashing_transcripts["Singlet", "transcripts"]
hashing@assays$RNA@data = as.matrix(log2(hashing@assays$RNA@counts + 1))

#Panel B Top, 1000 x 300
ribo_house = c(
  "RPL5",
  "RPL8",
  "RPL9",
  "RPL10A",
  "RPL11",
  "RPL14",
  "RPS5",
  "RPS6",
  "RPL13",
  "ACTB",
  "B2M",
  "GAPDH",
  "ALDOA",
  "PPIA"
)
DotPlot(hashing, features = ribo_house) +
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


genelist = rownames(hashing@assays$RNA@data)
ribolist = c(genelist[grep(c("^(?=.*RPL)(?!.*MRPL)"), genelist, perl = TRUE)], genelist[grep(c("^(?=.*RPS)(?!.*MRPS)"), genelist, perl =
                                                                                               TRUE)])
meanribo = colMeans(hashing@assays$RNA@data[ribolist, ], na.rm = TRUE)
hashing[["ribo"]] = meanribo

#Panel C Top, 400 x 300
vlncol = c("deepskyblue4", "deepskyblue3")
VlnPlot(
  hashing,
  features = c("ribo"),
  group.by = 'orig.ident',
  cols = vlncol,
  ncol = 1,
  pt.size = 0
) + ylab("log2 Expression") + xlab("") + NoLegend() + ggtitle("") +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 15,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#Fold-change calulation
hashing_ribo = data.frame(ribo = hashing@meta.data$ribo,
                          identity = hashing@meta.data$orig.ident)
hashing_ribo = aggregate(. ~ identity, hashing_ribo, median)
rownames(hashing_ribo) = hashing_ribo$identity
hashing_ribo["Doublet", "ribo"] / hashing_ribo["Singlet", "ribo"]


#Demuxlet PBMC Analysis
demuxlet = readRDS("demuxlet.subset")

#Panel A Right, 400 x 300
vlncol = c("hotpink4", "hotpink3")
VlnPlot(
  demuxlet,
  features = c("nCount_RNA"),
  group.by = 'status',
  cols = vlncol,
  ncol = 1,
  pt.size = 0
) + ylab("Transcripts") + xlab("") + NoLegend() + ggtitle("") +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 15,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#Fold-change calulation
demuxlet_transcripts = data.frame(
  transcripts = demuxlet@meta.data$nCount_RNA,
  identity = demuxlet@meta.data$status
)
demuxlet_transcripts = aggregate(. ~ identity, demuxlet_transcripts, median)
rownames(demuxlet_transcripts) = demuxlet_transcripts$identity
demuxlet_transcripts["Doublet", "transcripts"] / demuxlet_transcripts["Singlet", "transcripts"]

demuxlet@assays$RNA@data = as.matrix(log2(demuxlet@assays$RNA@counts + 1))
Idents(demuxlet) = "status"

#Panel B Bottom, 1000 x 300
ribo_house = c(
  "RPL5",
  "RPL8",
  "RPL9",
  "RPL10A",
  "RPL11",
  "RPL14",
  "RPS5",
  "RPS6",
  "RPL13",
  "ACTB",
  "B2M",
  "GAPDH",
  "ALDOA",
  "PPIA"
)
DotPlot(demuxlet, features = ribo_house) +
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
  ) # 1000 x 300

genelist = rownames(demuxlet@assays$RNA@data)
ribolist = c(genelist[grep(c("^(?=.*RPL)(?!.*MRPL)"), genelist, perl = TRUE)], genelist[grep(c("^(?=.*RPS)(?!.*MRPS)"), genelist, perl =
                                                                                               TRUE)]) #lookaround assertions to pull only RPL/RPS w/o mito genes
meanribo = colMeans(demuxlet@assays$RNA@data[ribolist, ], na.rm = TRUE)
demuxlet[["ribo"]] = meanribo

#Panel C Bottom, 400 x 300
vlncol = c("deepskyblue4", "deepskyblue3")
VlnPlot(
  demuxlet,
  features = c("ribo"),
  group.by = 'status',
  cols = vlncol,
  ncol = 1,
  pt.size = 0
) + ylab("log2 Expression") + xlab("") + NoLegend() + ggtitle("") +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 15,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#Fold-change calulation
demuxlet_ribo = data.frame(ribo = demuxlet@meta.data$ribo,
                           identity = demuxlet@meta.data$status)
demuxlet_ribo = aggregate(. ~ identity, demuxlet_ribo, median)
rownames(demuxlet_ribo) = demuxlet_ribo$identity
demuxlet_ribo["Doublet", "ribo"] / demuxlet_ribo["Singlet", "ribo"]



#12K Human-Mouse Mixture Analysis
mixture = readRDS("mixture.subset")

#Panel A Left, 400 x 300
vlncol = c("hotpink4", "hotpink3", "hotpink2") #pink!
VlnPlot(
  mixture,
  features = c("nCount_RNA"),
  group.by = 'orig.ident',
  cols = vlncol,
  ncol = 1,
  pt.size = 0,
  sort = TRUE
) + ylab("Transcripts") + xlab("") + NoLegend() + ggtitle("") +
  stat_summary(
    fun.y = median,
    geom = 'point',
    size = 15,
    colour = "grey27",
    shape = 95
  ) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#Fold-change calulation
mixture_transcripts = data.frame(
  transcripts = mixture@meta.data$nCount_RNA,
  identity = mixture@meta.data$orig.ident
)
mixture_transcripts = aggregate(. ~ identity, mixture_transcripts, median)
rownames(mixture_transcripts) = mixture_transcripts$identity
mixture_transcripts["Mixed", "transcripts"] / mixture_transcripts["Human", "transcripts"]
mixture_transcripts["Mixed", "transcripts"] / mixture_transcripts["Mouse", "transcripts"]
