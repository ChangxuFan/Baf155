rm(list = ls())
source("objects.R")
source("day9_only/step1.1a_helper_functions_2023-09-29.R")

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")

aoi$seurat_clusters <- sub("X", "", aoi$seurat_clusters)
plot.dir <- "day9_only/sth/step1.1"
plot.dir.rna <- "day9_only/sth/step1.1"

soi <- soi[, soi$seurat_clusters != "21"]

t.f.plot.umap(soi = soi, aoi = aoi, rna = T, atac = T, out.dir = "day9_only/sth/step1.1/umaps_rastered/",
              cell.type.map = "day9_only/data/cellType_map.tsv",
              color.map = "day9_only/data/colorMap.tsv")

markers <- c("Hlf", "Irf8", "Cebpe", "Epor", "Pf4", "Dntt")
plot.panel.list(panel.list = markers, obj = soi, order = F, assay = "RNA", label = F,
                root.name = paste0(plot.dir.rna, "/gex_markers_all.png"),
                publication = T, invisible = T, n.col = 2)
plot.panel.list(panel.list = markers, obj = soi, order = F, assay = "RNA", label = F,
                sample = SAMPLES.batch2,
                root.name = paste0(plot.dir.rna, "/gex_markers_d9.png"),
                publication = T, invisible = T, n.col = 2)


#>>>>>> BEGIN: bubble plot of markers for each cluster.

marker.df <- read.table("day9_only/data/cellType_markers.tsv", header = T,
                        sep = "\t")
m.bub <- marker.df$markers %>% strsplit(", ") %>% unlist()

lapply(c("h", "v"), function(x) {
  text.size = 10
  if (x == "h") {
    soi$cell_type <- factor(soi$cell_type, levels = rev(marker.df$cell_type))
    width <- 8
    height <- 4
  } else {
    soi$cell_type <- factor(soi$cell_type, levels = marker.df$cell_type)
    width <- 4
    height <- 8
  }
  
  if (x == "v") {
    m.bub <- rev(m.bub)
  }
  
  p <- DotPlot(soi, assay = "RNA", features = m.bub, cluster.idents = F, group.by = "cell_type", 
               dot.scale = 5, dot.min = 0.01) +
    theme() + xlab(NULL) + ylab(NULL)
  
  p <- p %>% utilsFanc::theme.fc.1(
    rotate.x.45 = T, italic.x = F, rm.x.ticks = F, text.size = text.size)
  
  if (x == "v") {
    p <- p + theme(axis.text.y = element_text(size = text.size, family = "Arial", face = "italic",
                                              color = "black"))
    p <- p + coord_flip()
  }
  
  p <- p + theme(legend.key.size = unit(0.2, "in"))
  ggsave(paste0("day9_only/sth/step1.1/bubblePlot_", x,".pdf"), p,  
         device = cairo_pdf, width = width, height = height, 
         dpi = 300)
  
})

