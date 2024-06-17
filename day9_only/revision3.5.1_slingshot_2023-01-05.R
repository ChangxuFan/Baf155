source("objects.R")
suppressMessages(library(slingshot))
wd <- "day9_only/sth/revision3.5.1_slingshot/"

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
sol <- lapply(SAMPLES, function(sample) {
  so <- soi[, soi$sample == sample]
})
names(sol) <- SAMPLES

slingos.repca <- utilsFanc::safelapply(SAMPLES, function(sample) {
  so <- sol[[sample]]
  so <- so[, so$seurat_clusters %in% as.character(c(0:6, 8:13))]
  # so <- so[, ! so$seurat_clusters %in% as.character(c(7, 16:21))]
  DefaultAssay(so) <- "RNA"
  assay <- "RNA"
  so <- NormalizeData(object = so, normalization.method = "LogNormalize", assay = assay)
  so <- FindVariableFeatures(so, assay = assay)
  so <- ScaleData(object = so, assay = assay)
  
  npc <- 30
  so <- RunPCA(object = so, assay = "RNA", npcs = npc)
  wd <- paste0(wd, "/rePCA_removeSmallClusters/")
  slingshot.seurat(so = so, cluster.ident = "seurat_clusters", start.clus = "1",
                   end.clus = c("10", "6", "12", "13", "3"),
                   approx_points = 400,
                   plot.dir = wd,
                   root.name = paste0("rePCA_removeSmallClusters_nPC", npc, "_", sample), 
                   saveRDS = T)
  
}, threads = 6)


slingo.list <- lapply(SAMPLES, function(sample) readRDS(paste0("day9_only/sth/revision3.5.1_slingshot/rePCA_removeSmallClusters/rePCA_removeSmallClusters_nPC30_", sample, "_slingshot.Rds")))
names(slingo.list) <- SAMPLES

lineage.list <- rep(list(1:5), 6)
names(lineage.list) <- SAMPLES

source("day9_only/step1.1a_helper_functions_2023-09-29.R")
soi <- soi[, soi$seurat_clusters != "21"]
soi <- t.f.cellType(
  soi, "day9_only/data/cellType_map.tsv")

color.df <- read.table("day9_only/data/colorMap.tsv",
                       comment.char = "", header = T, sep = "\t")
color.df <- color.df %>% filter(color != "")

cell.color.map <- soi@meta.data[, c("cell_type", "seurat_clusters"), drop = F] %>% 
  mutate(., cell.name = rownames(.)) %>% 
  left_join(color.df, by = "cell_type")

umap.list <- lapply(sol, function(so) so[, ! so$seurat_clusters %in% "21"]@reductions$umap@cell.embeddings)



sling.draw.curve(slingo.list = slingo.list, embed.to.list = umap.list, color.df = cell.color.map,
                 lineages.list = lineage.list, plot.dir = paste0(wd, "/rePCA_removeSmallClusters/"),
                 root = "drawCurve")

