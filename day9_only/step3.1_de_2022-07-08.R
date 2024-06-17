source("objects.R")
work.dir <- paste0("day9_only/sth/de_v1/")
dir.create(work.dir)
samples <- SAMPLES.batch2

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
soi <- soi[, soi$sample %in% samples]
soi$type[soi$type == "KO"] <- "cKO"

clusters <- as.character(c(0:6, 12:13))

de <- bulk.list(so = soi, assay = "RNA", slot = "counts", group.by = "sample",
                sample.order = samples, coldata.columns = c("type"), 
                threads = 1, design.formula = ~ type, contrast = c("type", "cKO", "WT"), 
                cluster.ident = "seurat_clusters", clusters = clusters, 
                filter.nz = T, filter.size = NULL, sequential.filter = T, 
                independentFiltering = T,
                single.sample = F, 
                do.plot = F, save.rds = T, 
                work.dir = work.dir, plot.dir = paste0(work.dir, "/plots/"))
de <- de[c(1,2,4)]

saveRDS(de, paste0(work.dir, "/bulk.list.small.Rds"))

de <- deseq2.summary(de, save = T, gene.out.dir = paste0(work.dir, "/summary/"))

genes.label <- list(c("Ccr2", "Itgam"), c("Esam", "Kit"), c("Gata2", "Epor"))
names(genes.label) <- names(de)

deseq2.xyplot(pbl = de, comp.vec = "cKO:WT", comp.meta = "type",
              transformation = log2, 
              plot.dir = "day9_only/sth/de_v1/plot_xy_label/", 
              device = "png", n.col = 3, theme.text.size = 15,
              add.label = T, label.list = genes.label, hjust = 0.5, vjust = 0.5, 
              use.repel = T, italic.label = T, y.lab = "KO", threads = 1)

#>>>>>>>> plot pathway enrichment for DEGs
msigdb.cats.short <- c("H..@C5..GO:BP")
de <- msigdb.m(de, summary.slot = "summary", 
               cats = msigdb.cats.short, threads = 4, universe = "expressed")

trash <- msigdb.sum(pbl = de, 
                    out.file = paste0("day9_only/sth/de_v1/msigdb_final/H_GOBP_expressed_summary.tsv"),
                    cats = msigdb.cats.short )
de <- de.change.name(de)
msigdb.dotplot(de, cat.to.plot = msigdb.cats.short, n = 10, 
               font.size = 10, 
               plot.out = paste0("day9_only/sth/de_v1/msigdb_final/H_GOBP_expressed_summary.png"))

