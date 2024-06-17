source("objects.R")
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")
wd <- "day9_only/sth/revision3.8.a_footPrint_pro_vs_enh/"
dir.create(wd)
motifPositions <- getPositions(ArchRProj = aoi, name = "homerFanc")
names(motifPositions)  <- make.unique(sub("/.+$", "", names(motifPositions)))

clusters <- c(0, 1, 3, 6, 12, 13) %>% as.character()

motifs <- list(
  "0" = c("IRF8.IRF", "^CEBP.bZIP.$"),
  "1" = "Etv2",
  "3" = c("Fli1", "Gata1"),
  "6" = c("IRF8.IRF", "^CEBP.bZIP.$"),
  "12" = c("E2A.bHLH.$"),
  "13" = c("Fli1", "Gata1")
)
greped <- lapply(motifs, function(motif) {
  names(motifPositions) %>% .[grepl(paste0(motif, collapse = "|"),.)]
}) 
names(greped) <- clusters

groups <- lapply(clusters, function(i) {
  groups <- aoi$cluster_type %>% unique() %>%
    .[grepl(paste0("^", i, "\\.\\."),.)] %>% sort()
  return(groups)
})
names(groups) <- clusters

#>>>>>>>> step2: separating by locations (promoter, intronic, intergenic): 
motifPositions <- lapply(motifPositions[unlist(greped)], function(gr)  {
  gr <- utilsFanc::gr.fast.annot(gr, "mm10", anno.cols = "annotation", use.klraps = F, 
                      simplify.location = T, narrow.promoter.boundaries = T, collapse.downstream = T)
})

location.levels <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic")
names(location.levels) <- location.levels

motifPositions.by.location <- lapply(motifPositions, function(gr) {
  split(gr, gr$annotation)
})
motifPositions.by.location <- lapply(location.levels, function(location) {
  grl <- lapply(motifPositions.by.location, function(x) {
    return(x[[location]])
  })
})
names(motifPositions.by.location) <- location.levels

aoi$cluster_type <- paste0(aoi$seurat_clusters, "..", aoi$type)
aoi <- addGroupCoverages(ArchRProj = aoi, groupBy = "cluster_type", useLabels = T)

# this will take a while to run:
seFoots.loc <- lapply(location.levels, function(location) {
  seFoots <- lapply(clusters, function(i) {
    seFoot <- getFootprints(ArchRProj = aoi, positions = motifPositions.by.location[[location]][greped[[i]]],
                            groupBy = "cluster_type", useGroups = groups[[i]])
  })
  names(seFoots) <- clusters
  return(seFoots)
})

saveRDS(seFoots.loc, paste0(wd, "/seFoots.loc.Rds"))

lapply(location.levels, function(location) {
  lapply(clusters, function(i) {
    lapply(c( "none"), function(norm) {
      footprint.decomposition(seFoot = seFoots.loc[[location]][[i]], out.dir = paste0(wd, "/byLocation/"),
                              root = paste0(location, "_", i, "_", norm), flankNorm = 50, debug = F)  
    })
  })
})

t.f.footprint.arrange(in.dir = paste0(wd, "/byLocation"), 
                      locations = c("Promoter", "Intron", "Intergenic"),
                      motifs.by.cluster = greped, plot.use = "foot",
                      out.dir = paste0(wd, "/byLocation_arrange_plot/"),
                      root = "foot")
