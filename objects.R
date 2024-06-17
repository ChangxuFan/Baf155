source("~/R_packages/common/base.R")
source("~/R_packages/common/sc.R")
suppressMessages(addArchRGenome("mm10"))
suppressMessages(library(DESeq2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GenomicInteractions))
suppressMessages(library(BiocParallel))
suppressMessages({
  R.files <- Sys.glob("R_*/*.R")
  for (R.file in R.files) {
    source(R.file)
  }
})

# all the rest of the .R files can be found at:
# https://github.com/ChangxuFan/common
# https://github.com/ChangxuFan/liteRnaSeqFanc
# https://github.com/ChangxuFan/scFanc

# these R scripts need to be sourced to use the functions I wrote to process the data.

SAMPLES <- c("WT_rep1", "KO_rep1", "WT_rep2", "KO_rep2", "WT_rep3", "KO_rep3")
# internal names. rep1: day 11 post 5FU; rep2/3: day 9 post 5FU

SAMPLES.batch2 <- c("WT_rep2", "KO_rep2", "WT_rep3", "KO_rep3")
SAMPLES.batch1 <- c("WT_rep1", "KO_rep1")

SAMPLES.hq <- c("WT_rep1", "KO_rep1", "WT_rep2", "WT_rep3", "KO_rep3")

SAMPLES.wt <- c("WT_rep1", "WT_rep2", "WT_rep3")
SAMPLES.ko <- c("KO_rep1", "KO_rep2", "KO_rep3")

SAMPLES.wt.batch2 <- c("WT_rep2", "WT_rep3")
SAMPLES.ko.batch2 <- c("KO_rep2", "KO_rep3")

EXON.ONLY.DIR <- "sth/count_exonOnly/" %>% normalizePath()
# results of cellranger count


# some metrics from cell ranger:
BC.METRICS.FILE.VEC.exonOnly <- paste0(EXON.ONLY.DIR, "/",  
                                       SAMPLES,
                                       "/outs/per_barcode_metrics.csv")

names(BC.METRICS.FILE.VEC.exonOnly) <- SAMPLES

BC.METRICS.FILE.LIST.exonOnly <- as.list(BC.METRICS.FILE.VEC.exonOnly)

FRAG.FILES.exonOnly <- paste0(EXON.ONLY.DIR, "/", SAMPLES, "/outs/atac_fragments.tsv.gz")

names(FRAG.FILES.exonOnly) <- SAMPLES

de.change.name <- function(de) {
  dict <- c("Neu", "HSPC", "Ery")
  names(dict) <- paste0("seurat_clusters_", c(0, 1,3))
  de <- lapply(de, function(s2b) {
    if (s2b$root.name %in% names(dict))
      s2b$root.name <- dict[s2b$root.name]
    return(s2b)
  })
  names(de) <- sapply(de, function(x) return(x$root.name))
  return(de)
}

