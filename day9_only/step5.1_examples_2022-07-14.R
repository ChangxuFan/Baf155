source("objects.R")
# This script generated the barplots of gene expression on the sides of the genome browser views
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
samples <- SAMPLES.batch2
soi <- soi[, soi$sample %in% samples]

# eg means example
eg.genes <- c("Itgam", "Ccr2", "Tnip3", "Esam", "Epor", "Car1", "Gata2", "Kit")

s2b.hm.single.gene(so = soi, group.ident.1 = "seurat_clusters", groups.1 = c("1", "0", "3"),
                   group.ident.2 = "sample", groups.2 = sort(samples),
                   assay = "RNA", slot = "data", 
                   genes = eg.genes, 
                   plot.dir = "day9_only/sth/de_v1/hm_ind/v2/")

