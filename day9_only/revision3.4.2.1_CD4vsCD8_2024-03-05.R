rm(list = ls())
source("objects.R")
wd <- "day9_only/sth/revision3.4.2.1_CD4vsCD8/"
dir.create(wd)
samples <- c("T.4.Nve.Sp.rpt1",
             "T.4.Nve.Sp.rpt2",
             "T.8.Nve.Sp.rpt1",
             "T.8.Nve.Sp.rpt2")
# these are from immgen

sample.info <- data.frame(sample = samples, 
                          celltype = paste0("CD", stringr::str_extract(samples, "4|8")),
                          peaks = paste0(wd, "/peaks/*", samples, "*narrowPeak"),
                          bam = paste0(wd, "/bam_open/*", samples, "*bam"),
                          project = "T")
sample.info$peaks <- Sys.glob(sample.info$peaks)
sample.info$bam <- Sys.glob(sample.info$bam)

# > sample.info
#            sample celltype                                                                    peaks
# 1 T.4.Nve.Sp.rpt1      CD4 day9_only/sth/revision3.4.2.1_CD4vsCD8//peaks/T.4.Nve.Sp.rpt1.narrowPeak
# 2 T.4.Nve.Sp.rpt2      CD4 day9_only/sth/revision3.4.2.1_CD4vsCD8//peaks/T.4.Nve.Sp.rpt2.narrowPeak
# 3 T.8.Nve.Sp.rpt1      CD8 day9_only/sth/revision3.4.2.1_CD4vsCD8//peaks/T.8.Nve.Sp.rpt1.narrowPeak
# 4 T.8.Nve.Sp.rpt2      CD8 day9_only/sth/revision3.4.2.1_CD4vsCD8//peaks/T.8.Nve.Sp.rpt2.narrowPeak
#                                                                         bam project
# 1 day9_only/sth/revision3.4.2.1_CD4vsCD8//bam_open/T.4.Nve.Sp.rpt1.open.bam       T
# 2 day9_only/sth/revision3.4.2.1_CD4vsCD8//bam_open/T.4.Nve.Sp.rpt2.open.bam       T
# 3 day9_only/sth/revision3.4.2.1_CD4vsCD8//bam_open/T.8.Nve.Sp.rpt1.open.bam       T
# 4 day9_only/sth/revision3.4.2.1_CD4vsCD8//bam_open/T.8.Nve.Sp.rpt2.open.bam       T

# these processed data can be easily generated using the AIAP pipeline.
# The .open.bam files were generated from the .open.bed files from the AIAP pipeline.
# They were converted to bam files using bedtools bedToBam

peaks <- peak.merge(peak.files = sample.info$peaks, score.col = 9, th = 8, threads = 6)

peakmat <- bam.count(sample.info = sample.info, bam.col = "bam", sample.col = "sample",
                     features.gr = peaks, bSingleEnd = T, sort = T, return.mat = T, threads = 6)
saveRDS(peaks, paste0(wd, "/peakset.Rds"))
saveRDS(peakmat, paste0(wd, "/peakmat.Rds"))
write.table(sample.info, paste0(wd, "/sample_info.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)

da <- diffbind.2.raw.a2bl(sample.info = sample.info, mat = peakmat,
                          cluster.ident = "project",sample.col = "sample", 
                          filter.nz = F, filter.size = NULL, sequential.filter = T,
                          coldata.columns = c("celltype"), 
                          design.formula = ~ celltype,
                          contrast = c("celltype", "CD8", "CD4"),
                          threads = 1,
                          work.dir = wd)

da <- deseq2.summary(da, gene.out.dir = paste0(wd, "/summary/"))

deseq2.xyplot(pbl = da, comp.vec = "CD8:CD4", comp.meta = "celltype", 
              transformation = function(x) log2(x + 1), 
              plot.dir = paste0(wd, "/plot_xy/"), root.name = "log2",
              device = "png", n.col = length(da), theme.text.size = 15, 
              plot.each.comp = F)

dasc <- readRDS("day9_only/sth/da_v1/ao.bulk.list.Rds")

lympho.peaks <- rownames(dasc$seurat_clusters_12$bulk.mat) %>% utilsFanc::loci.2.gr()

#>>>>>>>
# The following code was used to generate the box plot of accessibility fold changes.
# The function da.motif.foldChange() was originally written to compare the fold changes
# of peaks with different motifs. Here, it was adopted to compare CD4 vs CD8 peaks, hence "pseudo"

pseudo.peak.motif.df <- lapply(c("up", "down"), function(type) {
  dar <- da[["T"]]$summary[[paste0(type, ".genes")]] %>% utilsFanc::loci.2.gr()
  dar$CD48_loci <- utilsFanc::gr.get.loci(dar)
  dar$motif <- ifelse(type == "up", "CD8", "CD4")
  dar$score <- 10
  dar$n.motif <- 1
  
  lympho.peaks$gene <- utilsFanc::gr.get.loci(lympho.peaks)
  dar <- plyranges::join_overlap_left(dar, lympho.peaks)
  dar <- dar[!is.na(dar$gene)]
  df <- mcols(dar) %>% as.data.frame()
}) %>% do.call(rbind, .)

write.table(pseudo.peak.motif.df, paste0(wd, "/pseudo.peak.motif.df.tsv"), 
            row.names = F, col.names = T, quote = F, sep = "\t")

df <- da.motif.foldChange(da = dasc["seurat_clusters_12"], comps = 'beta', peak.motif.df = pseudo.peak.motif.df,
                    plot.dir = paste0(wd, "/peak_fc/"), motifs.use = c("CD8", "CD4"),
                    root.name = "CD8vsCD4", bViolin = F, bBoxplot = T, bDot = F, write.each.df = T,
                    publication = T, pt.size = 0.5, alpha = 0.8)
