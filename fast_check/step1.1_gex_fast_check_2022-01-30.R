source("objects.R")
work.dir <- "fast_check/sth/rna/"
plot.dir <- "fast_check/plots/"

paste0("mkdir -p ", plot.dir, " ", work.dir) %>% system()

samples <- SAMPLES

sample.info <- data.frame(
  sample = samples,
  dir = EXON.ONLY.DIR %>% paste0("/", samples, "/outs/filtered_feature_bc_matrix"),
  type = samples %>% sub("_rep.+", "",.),
  batch = c("batch1", "batch1", rep("batch2", 4)),
  rep = samples %>% sub(".+_", "", .)
)

write.table(x = sample.info, "fast_check/sample_info.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)

sol <- sc.qc.construct.so(sample.info = "fast_check/sample_info.tsv",
                          project.name = "fast_check",
                          mt.pattern = "mt-")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      override = T)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "ori")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = rep(list(750), 6) %>% `names<-`(samples),
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end")

sol <- sc.qc.elbow.filter(sol = sol,metas = c("nFeature_RNA"),
                          project.name = "fast_check", take.lower = F)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end_post")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = list(WT_rep1 = 3700, KO_rep1 = 3700,
                                                        WT_rep2 = 3000, KO_rep2 = 3700,
                                                        WT_rep3 = 3000, KO_rep3 = 3700),
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "higher_end")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "percent.mt",
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "percent.mt", ymax = 100, plot.dir = plot.dir)

sol <- sc.qc.elbow.filter(sol = sol,metas = c("percent.mt", "nFeature_RNA"),
                          project.name = "fast_check")

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000, plot.dir = plot.dir, project.name = "post_filter")
t <- sc.rank.qc(sol = sol, feature = "percent.mt", ymax = 100, plot.dir = plot.dir, project.name = "post_filter")


sol <- utilsFanc::safelapply(sol, function(so) {
  root.name <- so@meta.data$sample[1]
  utilsFanc::t.stat(root.name)
  so <- cluster.pipe(soi = so, assay = "SCT", pc.dim = 1:30, cluster.resolution = 0.7, 
                     work.dir = paste0(work.dir, "/per_sample/", root.name),
                     plot.dir = paste0(work.dir, "/per_sample/", root.name, "/plots/"), 
                     project.name = root.name, metas.include = "sample", 
                     sample.tsv = "fast_check/sample_info.tsv",
                     save.rds = T, do.sct = T, sct.n.vars = 3000, hm.dims = 11:30, 
                     plot.common.markers = T, bc.metrics.file.list = BC.METRICS.FILE.LIST.exonOnly)
  return(so)
}, threads = length(samples))
