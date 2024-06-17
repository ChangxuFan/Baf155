source("objects.R")
wd <- "day9_only/sth/revision3.8.d_more_clusters"
dir.create(wd)

de <- readRDS("day9_only/sth/de_v1/bulk.list.Rds")
de <- de[paste0("seurat_clusters_", c(6, 12, 13))]
de <- deseq2.summary(de, gene.out.dir = paste0(wd, "/de_summary"))

deseq2.xyplot(pbl = de, comp.vec = "cKO:WT", comp.meta = "type",
              transformation = function(x)log2(x + 1), 
              plot.dir = wd, root.name = "de",
              device = "png", n.col = 3, theme.text.size = 15,
              add.label = F, 
              y.lab = "KO", threads = 1)

da <- readRDS("day9_only/sth/da_v1/ao.bulk.list.Rds")
da <- da[paste0("seurat_clusters_", c(6, 12, 13))]
da <- deseq2.summary(da, gene.out.dir = paste0(wd, "/da_summary"))

deseq2.xyplot(pbl = da, comp.vec = "cKO:WT", comp.meta = "type",
              transformation = function(x)log2(x + 1), 
              plot.dir = wd, root.name = "da",
              device = "png", n.col = 3, theme.text.size = 15,
              add.label = F, 
              y.lab = "KO", threads = 1)

#>>>>>>>> homer:
samples.bg <- SAMPLES.wt.batch2
scatter.df <- data.frame(x = SAMPLES.wt.batch2, y = SAMPLES.ko.batch2)

a2bl.homer.pipe(a2bl = da, work.dir = paste0(wd, "/homer"),
                slot = "res.exp", rank.key = "padj", 
                desc = F, top.n.vec = c(1000), n.bg.each = 25, no.replace = T,
                samples.for.bg = samples.bg, 
                scatter.xy.df = scatter.df, 
                threads.a2bl = 1, threads.titrate = 1, threads.homer = 10, 
                genome = "mm10", size = "given", denovo = F, homer.run = T)

#>>>>>>>> Print DEG and DAR
da <- readRDS("day9_only/sth/da_v1/ao.bulk.list.Rds")
da <- da[paste0("seurat_clusters_", c(0, 1, 3, 6, 12, 13))]
da <- deseq2.summary(da)
alias <- data.frame(name = paste0("seurat_clusters_", c(0, 1, 3, 6, 12, 13)),
                    alias = c("Neutrophil progenitors", "HSPC", "Erythroid progenitors", "Mono_Mac_DC progenitors",
                              "Lymphoid progenitors", "Megakaryocytic progenitors"))
de.write.xls(da, cluster.alias = alias, out.file = paste0(wd, "/DARs.xlsx"))


de <- readRDS("day9_only/sth/de_v1/bulk.list.Rds")
de <- de[paste0("seurat_clusters_", c(0, 1, 3, 6, 12, 13))]
de <- deseq2.summary(de)
alias <- data.frame(name = paste0("seurat_clusters_", c(0, 1, 3, 6, 12, 13)),
                    alias = c("Neutrophil progenitors", "HSPC", "Erythroid progenitors", "Mono_Mac_DC progenitors",
                              "Lymphoid progenitors", "Megakaryocytic progenitors"))
de.write.xls(de, cluster.alias = alias, out.file = paste0(wd, "/DEGs.xlsx"))
