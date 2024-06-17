source("objects.R")
work.dir <- "sync_all/sth/atac/merge_v1/"
plot.dir <- paste0(work.dir, "/Plots/")
system(paste0("mkdir -p ", work.dir))

samples <- SAMPLES

aoi <- readRDS(paste0(work.dir, "/ao_IterativeLSI.Rds"))
aoi <- archr.macs2.pipe(ao = aoi, work.dir = work.dir, cluster.ident = "seurat_clusters", 
                        peak.annot = "homer")
