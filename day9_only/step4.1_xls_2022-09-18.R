source("objects.R")
de <- readRDS("day9_only/sth/de_v1/bulk.list.small.Rds")
de <- deseq2.summary(de)
alias.df <- data.frame(name = paste0("seurat_clusters_", c(0, 1, 3)),
                       alias = c("Neu", "HSPC", "Ery"))
de.write.xls(de = de[1:3], cluster.alias = alias.df,
             out.file = "day9_only/xlsx/DEG_day9.xlsx")
da <- readRDS("day9_only/sth/da_v1/ao.bulk.list.small.Rds")
da <- deseq2.summary(da)
de.write.xls(de = da[1:3], cluster.alias = alias.df,
             out.file = "day9_only/xlsx/DAR_day9.xlsx")
