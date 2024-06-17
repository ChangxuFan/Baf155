# This script shows that peaks with IRF8 motif is less affected by Baf155 KO, compared to
# those with the CEBP motif.
source("objects.R")
wd <- "day9_only/sth/revision3.7.1_mac_vs_neu/"
dir.create(wd)
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")

homer.motifs <- abaFanc2::homer.motif.2.pwMatrix("~/motifs/homer/sth/homer_motif_database.txt")
# You can download this from http://homer.ucsd.edu/homer/custom.motifs
aoi <- addMotifAnnotations(ArchRProj = aoi, name = "homerFanc", motifPWMs = homer.motifs)
saveRDS(aoi, "sync_all/sth/atac/merge_v1/ao.Rds")

da <- readRDS("day9_only/sth/da_v1/ao.bulk.list.Rds")

clusters <- paste0("seurat_clusters_", c(0, 1, 2, 3, 6))

motifPositions <- getPositions(ArchRProj = aoi, name = "homerFanc")
motifs <- c("CEBP.bZIP./ThioMac", "IRF8.IRF./")
greped <- names(motifPositions) %>% .[grepl(paste0(motifs, collapse = "|"),.)]

peak.motif.df <- archr.peak.motif.df.gen(motifPositions[greped], aoi = aoi, threads = 1, 
                                         out.rds = paste0(wd, "/peak_motif_df_CEBP_IRF8.Rds"))

da.motif.foldChange(da = da[clusters], comps = 'beta', peak.motif.df = peak.motif.df,
                    plot.dir = paste0(wd, "/motif_fc/"), motifs.use = greped,
                    root.name = "CEBP_IRF8", bViolin = F, bBoxplot = T, bDot = F, write.each.df = F,
                    publication = T)

#<<<<<<<<END: boxplot of log2FC

#>>>>>>>>BEGIN: plot homer results
motifs <- c("CEBP(bZIP)/ThioMac-CEBPb-ChIP-Seq(GSE21512)/Homer", "IRF8(IRF)/BMDM-IRF8-ChIP-Seq(GSE77884)/Homer")
names(motifs) <- c("CEBP", "IRF8")
homer.txt <- "day9_only/sth/revision3.8.d_more_clusters/homer/seurat_clusters_6/res.exp_padj_top_1000_down/homer/knownResults.txt"
homer.plot.pct(homer.txt = homer.txt, motifs = motifs, rm.legend = T, 
               out.dir = paste0(wd, "/motif_pct_homer/"))
