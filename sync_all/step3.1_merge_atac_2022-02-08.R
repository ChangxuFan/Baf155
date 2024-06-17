source("objects.R")
work.dir <- "sync_all/sth/atac/merge_v1/"
plot.dir <- paste0(work.dir, "/Plots/")
system(paste0("mkdir -p ", work.dir))

samples <- SAMPLES

arrow.files <- paste0("fast_check/sth/atac//", samples, "/", samples, ".arrow")


aoi <- ArchRProject(ArrowFiles = arrow.files, outputDirectory = work.dir,
                    copyArrows = F)
# To the reader: you might want to use copyArrows = T, since the arrow files provided are linked to my directory structure...

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")

cells <- intersect(aoi$cellNames, colnames(soi))
soi <- soi[, cells]
aoi <- aoi[cells, ]

saveRDS(aoi, paste0(work.dir, "/ao.Rds"))


trash <- archr.cluster.pipe(ao = aoi, minTSS = 10, bc.metrics.file.list = NULL,
                            dimsToUse.umap = 1:15, dimsToUse.cluster = 1:10, cluster.res = 0.6,
                          work.dir = work.dir, force = T, do.harmony = F, sample.order = samples)

