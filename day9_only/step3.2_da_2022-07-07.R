source("objects.R")

aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")
aoi$type[aoi$type == "KO"] <- "cKO"
aoi$seurat_clusters <- sub("X", "", aoi$seurat_clusters)

samples <- SAMPLES.batch2
filter.samples <- intersect(SAMPLES.batch2, SAMPLES.hq)
samples.bg <- SAMPLES.wt.batch2

scatter.df <- data.frame(x = SAMPLES.wt.batch2, y = SAMPLES.ko.batch2)
work.dir <- "day9_only/sth/da_v1/"
dir.create(work.dir, showWarnings = F, recursive = T)

peakmat <- readRDS("sync_all/sth/atac/merge_v1/peakmat_full.Rds")
peakmat <- assay(peakmat)


cluster.ident <- "seurat_clusters"
clusters <- as.character(c(0:6, 12:13))
filter.size <- 6

da <- ao.2.bulk.list(ao = aoi, peakmat = peakmat, pseudo = F,
                     filter.samples = filter.samples,
                     filter.nz = T, filter.size = c(filter.size, 0.98), sequential.filter = T,
                     independentFiltering = T,
                     cell.prep.dir = paste0(work.dir, "/cell_prep/"),
                     coldata.columns = c("type"),
                     cluster.ident = cluster.ident, clusters = clusters,
                     group.ident = "Sample", groups = samples,
                     threads = 1, design.formula = ~ type, contrast = c("type", "cKO", "WT"),
                     work.dir = work.dir, do.plot = F)
da <- readRDS(paste0(work.dir, "/ao.bulk.list.Rds"))
da <- deseq2.summary(pbl = da, padj.cutoff = 0.05, log2fc.cutoff = 1,
                     gene.out.dir = paste0(work.dir, "/summary/"),
                     stats.file = paste0(work.dir, "/summary/stats.tsv"))

da <- da[c(1,2,4,5)]
saveRDS(da, paste0(work.dir, "/ao.bulk.list.small.Rds"))
# 
deseq2.xyplot(pbl = da, comp.vec = "cKO:WT", comp.meta = "type",
              transformation = log2,
              plot.dir = paste0(work.dir, "/plot_xy/"), root.name = "log2",
              device = "png", n.col = length(da), theme.text.size = 15)

#>>>>>>>> motif enrichment in DARs:

a2bl.homer.pipe(a2bl = da, work.dir = paste0(work.dir, "/homer"),
                slot = "res.exp", rank.key = "padj", 
                desc = F, top.n.vec = c(1000), n.bg.each = 25, no.replace = T,
                samples.for.bg = samples.bg, 
                scatter.xy.df = scatter.df, 
                threads.a2bl = 4, threads.titrate = 2, threads.homer = 1, 
                genome = "mm10", size = "given", denovo = F, homer.run = T)

