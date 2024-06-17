source("objects.R")
#>>>>>>>> BEGIN: plot chromVAR
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")

samples <- SAMPLES.batch2
aoi <- aoi[aoi$Sample %in% samples, ]

aoi$seurat_clusters_ori <- aoi$seurat_clusters %>% sub("X", "", .) %>% 
  factor(., levels = sort(unique(as.numeric(.))))

# Only keeping large clusters
aoi.rm <- aoi[aoi$seurat_clusters_ori %in% as.character(0:13),]

# Remove the cycling cluster:
aoi.no7 <- aoi.rm[aoi.rm$seurat_clusters_ori != "7", ]


ao <- aoi.no7
ao <- addImputeWeights(ao)

motifs <- c("GATA1", "Klf1" ,"CEBP", "EBF1", "IRF4", "TBET", "PAX5", "FLI1", "RUNX1")
markerMotifs <- getFeatures(ao, select = paste(motifs, collapse="|"), useMatrix = "homerMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

p <- plotEmbedding(
  ArchRProj = ao, 
  colorBy = "homerMatrix", 
  name = markerMotifs, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(ao)
)

trash <- wrap.plots.fanc(plot.list = p, plot.out = "day9_only//chromVAR/motifs_v1.pdf")
#>>>>>>>> END: plot chromVAR

#>>>>>>>> BEGIN: plot marker motifs
markerPeaks <- readRDS("sync_all/sth/atac/merge_v1/markerPeaks_seurat_clusters.Rds")
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = ao,
  peakAnnotation = "homer",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

clusters.use <- paste0("X", c(0, 1, 3, 6, 12, 13))

motifs.selected <- c("Gata1", "EKLF", "CEBP.bZIP_34", "IRF4", "PU.1.ETS", "RUNX1", "E2A.bHLH_51",
                     "Tcf4", "HOXA9", "FLI1", "Mef2c", "STAT1")

motifs.selected.se <- enrichMotifs %>% .[grep(paste0(motifs.selected, collapse = "|"),
                                              rownames(.)), ]
rownames(motifs.selected.se)[grepl("FLI1", rownames(motifs.selected.se))] <- "FLI1"
rownames(motifs.selected.se) <- rownames(motifs.selected.se) %>% sub("\\..+$", "", .)

save.base.plot(exp = plotEnrichHeatmap(motifs.selected.se[, clusters.use],
                                       n = 30, transpose = F, cutOff = 10),
               file = "day9_only/motif_enrich/homer_selected.pdf",
               height = 300, width = 600)

##########
# notes on the TFs selected: They are supported by literature
# Mef2c: https://pubmed.ncbi.nlm.nih.gov/19211936/
# Stat1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2096448/
# Hoxa9: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1895111/

#<<<<<<<< BEGIN: plot marker motifs
