source("objects.R") 
work.dir <- "sync_all/sth/rna/int_rpca20"
plot.dir <- paste0(work.dir, "/plots/")
dir.create(plot.dir, showWarnings = F, recursive = T)

samples <- SAMPLES

sol <- utilsFanc::safelapply(samples, function(sample) {
  so <- readRDS(paste0("sync_all/sth/rna","/", sample, "/soi.Rds"))
  return(so)
}, threads = length(samples))

names(sol) <- samples
features <- SelectIntegrationFeatures(object.list = sol, nfeatures = 3000)
sol <- PrepSCTIntegration(object.list = sol, anchor.features = features)
sol <- lapply(X = sol, FUN = RunPCA, features = features)
names(sol) <- samples

anchors <- FindIntegrationAnchors(object.list = sol, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, 
                                         reduction = "rpca", k.anchor = 20)
soi.rpca <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
saveRDS(soi.rpca, paste0(work.dir, "/soi_v1.Rds"))

soi.rpca <- cluster.pipe(
  soi = soi.rpca, assay = "integrated",
  pc.dim = 1:30, cluster.resolution = 0.8,
  ####
  #### perhaps should change the resolution to 0.8
  ####
  work.dir = work.dir, plot.dir = plot.dir,
  project.name = "rpca", metas.include = "sample",
  save.rds = T, do.sct = F,
  plot.common.markers = T, split.by = "sample", split.order = samples
)
# soi.rpca <- readRDS(paste0(work.dir, "/soi.Rds"))
int.corr(q.vec = sol, s = soi.rpca, meta = "seurat_clusters", out.file = paste0(work.dir, "/int_cor.tsv"))
int.corr.plot(soi = soi.rpca, sol = sol, cluster.ident = "seurat_clusters", order = F,
              plot.dir = paste0(work.dir, "/int_corr_plot/"))
