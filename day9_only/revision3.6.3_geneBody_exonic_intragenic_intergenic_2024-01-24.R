source("objects.R")
da <- readRDS("day9_only/sth/da_v1/ao.bulk.list.small.Rds")

a2bl.distro.titrate.log2FC(a2bl = da[1:3], plot.out = "day9_only/sth/revision3.6.3_proEnh/proEnh/proEnh.png",
                           width = 6, text.size = 14, cluster.levels = paste0("seurat_clusters_", c(3, 1, 0)))


a2bl.distro.titrate.log2FC(a2bl = da[1:3], mode = "location", ### key change
                           log2FC.vec = c(0, 0.5, 1, 1.5, 2),
                           plot.out = "day9_only/sth/revision3.6.3_proEnh/exonIntronInter/exonIntronInter.png",
                           width = 6, text.size = 14, cluster.levels = paste0("seurat_clusters_", c(3, 1, 0)))

