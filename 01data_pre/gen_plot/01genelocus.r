library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)

#========================================================
# <U+1F525>关键1：关闭 UCSC 限制（必须最前）
#========================================================

options(ucscChromosomeNames = FALSE)

#========================================================
# 参数
#========================================================

chr0 <- "3"
pos  <- 95196819
flank <- 250000

from <- pos - flank
to   <- pos + flank

gtf_file <- "../Sus_scrofa.Sscrofa11.1.100.gtf"

#========================================================
# 读取 GTF
#========================================================

gtf <- import(gtf_file)

#========================================================
# <U+1F525>关键2：过滤 scaffold（必须）
#========================================================

gtf <- gtf[grepl("^[0-9XY]+$", seqnames(gtf))]

#========================================================
# 只保留 gene
#========================================================

genes <- gtf[gtf$type == "gene"]

#========================================================
# 区域筛选
#========================================================

genes_sub <- genes[
  seqnames(genes) == chr0 &
    start(genes) <= to &
     end(genes) >= from
]

cat("Genes found:", length(genes_sub), "\n")

if (length(genes_sub) == 0) {
  stop("<U+274C> No genes found in region.")
}

#========================================================
# gene name
#========================================================

genes_sub$symbol <- genes_sub$gene_name
genes_sub$symbol[is.na(genes_sub$symbol)] <- genes_sub$gene_id[is.na(genes_sub$symbol)]

#========================================================
# GeneRegionTrack
#========================================================

grtrack <- GeneRegionTrack(
  range = genes_sub,
  genome = "susScr11",
  chromosome = chr0,
  name = "Genes",

  transcriptAnnotation = "symbol",
  geneSymbol = TRUE,
  showId = TRUE,

  collapseTranscripts = "meta",
  height = 0.3,        # <U+2B50>控制“整条轨道高度”（最重要）
  stackHeight = 0.3,   # <U+2B50>控制每条transcript高度
  min.height = 0.2,    # <U+2B50>防止太粗

  fill = "grey70",
  col = "black"
)

#========================================================
# SNP track
#========================================================

snpTrack <- AnnotationTrack(
  start = pos,
  end = pos,
  chromosome = chr0,
  genome = "susScr11",
  name = "Lead SNP",
  fill = "red",
  col = "red"
)

axisTrack <- GenomeAxisTrack()

#========================================================
# 输出
#========================================================

pdf("Gviz_fixed_final.pdf", width = 10, height = 4)

plotTracks(
  list(axisTrack, grtrack, snpTrack),
  from = from,
  to = to,
  background.panel = "white"
)

dev.off()

cat("Done\n")
