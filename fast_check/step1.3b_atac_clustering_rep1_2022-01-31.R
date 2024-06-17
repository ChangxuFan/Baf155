source("objects.R")
# rep1 (day 11) samples were processed first (this script). They were processed slightly differently. 
# after obtainind day 9 samples (rep2/3), we re-ran most of the analyses for rep1 together with them,
# with the exception of doublet removal.
# This script is added so that the subsequent doublet removal scripts would make sense.

work.dir <- "fast_check/sth/rep1/rna/"
plot.dir <- paste0(work.dir, "/plots/")

paste0("mkdir -p ", work.dir, " ", plot.dir) %>% system()

samples <- c("WT_rep1", "KO_rep1")

sample.info <- read.table("fast_check/sample_info.tsv", header = T) %>% 
  filter(sample %in% samples)

sol <- sc.qc.construct.so(sample.info = sample.info,
                          project.name = "fast_check",
                          mt.pattern = "mt-")


sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      override = T)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "ori")

# need to filter out the lower end
sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = list(WT_rep1 = 750, KO_rep1 = 750),
                                      override = T)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end")

sol <- sc.qc.elbow.filter(sol = sol,metas = c("nFeature_RNA"),
                          project.name = "fast_check", take.lower = F)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end_post")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = list(WT_rep1 = 3000, KO_rep1 = 3500),
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

saveRDS(sol, paste0(work.dir, "/sol_v1.Rds"), compress = F)

# sol <- readRDS( paste0(work.dir, "/sol_v1.Rds"))

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
}, threads = 2)



addArchRThreads(threads = 4) 
addArchRThreads(threads = 12) 
work.dir <- "fast_check/sth/rep1/"

system("mkdir -p " %>% paste0(work.dir))
samples <- c("WT_rep1", "KO_rep1")

frag.files <- FRAG.FILES.exonOnly[samples]

aol <- utilsFanc::safelapply(names(frag.files), function(sample) {
  frag.file <- frag.files[sample]
  sub.work.dir <- paste0(work.dir, "/", sample)
  ao <- ao.gen(frag.files = frag.file, minTSS = 10, minFrags = 1000,
               work.dir = sub.work.dir)
  saveRDS(ao, paste0(sub.work.dir, "/ao.Rds"))
  return(ao)
}, threads = 1)



samples <- c("WT_rep1", "KO_rep1")
addArchRThreads(threads = 12) 
work.dir <- "fast_check/sth/rep1/"

work.dir.rna <- "fast_check/sth/rep1/rna/per_sample/"
system("mkdir -p " %>% paste0(work.dir))

utilsFanc::safelapply(samples, function(sample) {
  work.dir <- paste0(work.dir, "/", sample)
  work.dir.rna <- paste0(work.dir.rna, "/", sample)
  ao <- readRDS(paste0(work.dir, "/ao.Rds")) %>% recoverArchRProject()
  so <- readRDS(paste0(work.dir.rna, "/soi.Rds"))
  ao <- ao[ao$cellNames %in% get.cell.names.seurat(so = so, style = "ArchR"),]
  ao <- archr.cluster.pipe(ao = ao, so = so, minTSS = 10, work.dir = work.dir, do.harmony = F,
                           force = T)
  return()
}, threads = 1)

work.dir <- "fast_check/sth/rep1/"
samples <- c("WT_rep1", "KO_rep1")

aol.fu <- lapply(samples, function(x) {
  readRDS(paste0(work.dir, "/", x, "/ao_IterativeLSI.Rds")) %>% 
    return()
})
names(aol.fu) <- samples

lapply(samples, function(x) {
  archr.titrate.doublets(ao = aol.fu[[x]],
                         plot.dir = paste0(work.dir, "/", x, "/doublets/"),
                         root.name = "ao", 
                         label.vars = c("Clusters", "seurat_clusters"),
                         remove.frac = c(0.04, 0.08, 0.10))
})


semi.clean.rmClus.list <- list(
  WT_rep1 = paste0("C", c(2, 4, 11)),
  KO_rep1 = paste0("C", c(1, 3, 10, 11))
)

semi.clean.rmSeuratClus.list <- list(
  WT_rep1 = c("13", "14"),
  KO_rep1 = c("5", "7", "11")
)

round.2.pct <- list(
  WT_rep1 = c(0.02, 0.03, 0.04),
  KO_rep1= c(0.02, 0.03, 0.04)
)

stats.df <- utilsFanc::safelapply(samples, function(sample) {
  ao <- aol.fu[[sample]]
  semi.dir <- paste0(work.dir, "/", sample, "/semiClean/")
  dir.create(semi.dir, showWarnings = F, recursive = T)
  
  stats.df <- data.frame(sample = sample, pre.remove = ao$cellNames %>% length())
  
  ao <- ao[! ao$Clusters %in% semi.clean.rmClus.list[[sample]],]
  ao <- ao[! ao$seurat_clusters %in% semi.clean.rmSeuratClus.list[[sample]], ]
  
  saveRDS(ao, paste0(semi.dir, "/ao.Rds"))
  stats.df$post.remove <- ao$cellNames %>% length()
  write.table(stats.df, paste0(semi.dir, "/stats.tsv"), quote = F, col.names = T,
              row.names = F, sep = "\t")
  archr.titrate.doublets(ao = ao, plot.dir = semi.dir, root.name = "ao_semiClean",
                         label.vars = c("Clusters", "seurat_clusters"),
                         remove.frac = round.2.pct[[sample]])
  return(stats.df)
}, threads = 1) %>% Reduce(rbind, .)

write.table(stats.df, paste0(work.dir, "/fu_semiClean_stats.tsv"), quote = F, col.names = T,
            row.names = F, sep = "\t")

round.2.pct <- c(0.04, 0.04)
names(round.2.pct) <- samples

utilsFanc::safelapply(samples, function(sample) {
  cleaned.v1.dir <- paste0(work.dir, "/", sample, "/cleaned_v1/")
  dir.create(cleaned.v1.dir, showWarnings = F, recursive = T)
  ao.source <- paste0(work.dir,"/", sample, "/semiClean/ao_semiClean_",round.2.pct[sample], "_ao.Rds")
  cmd <- paste0("cp ", ao.source, " ", cleaned.v1.dir, "/ao.Rds")
  
  print(cmd); system(cmd)
}, threads = 1)
