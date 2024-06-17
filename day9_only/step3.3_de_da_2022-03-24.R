source("objects.R")
# this script is to show that downDEGs usually have reduced promoter chromatin accessibility, but upDEGs
# usually do not have increased promoter chromatin accessibility.

de <- readRDS("day9_only/sth/de_v1/bulk.list.small.Rds")
de <- deseq2.summary(de)

da <- readRDS("day9_only/sth/da_v1/ao.bulk.list.small.Rds")
da <- deseq2.summary(da)

# gtf <- readRDS("~/genomes/mm10/gencode/gencode.vM25.annotation.gff3.Rds")
gtf <- rtracklayer::import("~/genomes/mm10/gencode/gencode.vM25.annotation.gff3.gz")
utilsFanc::safelapply(c(1), function(fc) {
  p <- de.enrich.da(de[c(1:3)], da, gtf, da.fc.filter = fc,
                    plot.out = paste0("day9_only/de_da/de_da_fc_", fc, ".png"))
  return()
}, threads = 1)

# let's make sure there is nothing wrong with the code.

de.enrich.da.check(de = de[c(1:3)], da = da,
                   tss.file = "~/genomes/mm10/TSS/gencode.vM25.TSS.bed", 
                   out.dir = "day9_only/de_da/doubleCheck/")
# I made sure they match. Also checked a few examples on scViewer2 and genome browser:
# Bmyc Gramd1b Olr1 Olfr417 