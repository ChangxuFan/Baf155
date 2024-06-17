source("objects.R")
## change addArchRthreads if you change the number of samples to parallel
addArchRThreads(threads = 3) 
work.dir <- "fast_check/sth/atac/"

work.dir.rna <- "fast_check/sth/rna/per_sample/"
system("mkdir -p " %>% paste0(work.dir))

samples <- SAMPLES

utilsFanc::safelapply(samples, function(sample) {
  work.dir <- paste0(work.dir, "/", sample)
  work.dir.rna <- paste0(work.dir.rna, "/", sample)
  ao <- readRDS(paste0(work.dir, "/ao.Rds")) %>% recoverArchRProject()
  so <- readRDS(paste0(work.dir.rna, "/soi.Rds"))
  ao <- ao[ao$cellNames %in% get.cell.names.seurat(so = so, style = "ArchR"),]
  ao <- archr.cluster.pipe(ao = ao, so = so, minTSS = 10, work.dir = work.dir, do.harmony = F,
                           force = T, dimsToUse = 1:15, dimsLSI = 1:15)
  return()
}, threads = 6)


