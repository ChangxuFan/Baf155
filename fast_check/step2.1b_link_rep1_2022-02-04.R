source("objects.R")
batch1.dir <- "fast_check/sth/rep1/"
batch2.dir <- "fast_check/sth/atac/"
lapply(SAMPLES.batch1, function(sample) {
  ao.old <- paste0(batch1.dir, "/", sample, "/cleaned_v1/ao.Rds") %>% readRDS()
  ao.new <- paste0(batch2.dir, "/", sample, "/ao_IterativeLSI.Rds") %>% readRDS()
  cells.old <- ao.old$cellNames
  
  ao.new <- ao.new[ao.new$cellNames %in% cells.old, ]
  
  out.dir <- paste0(batch2.dir, "/", sample ,"/cleaned_v1/")
  dir.create(out.dir, showWarnings = F)
  saveRDS(ao.new, paste0(out.dir, "/ao.Rds"))
  write(length(ao.new$cellNames), paste0(out.dir, "/cell_number.txt"))
  return()
})