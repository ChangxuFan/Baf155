source("objects.R")
wd <- "day9_only/sth/revision1.4_markers/"
dir.create(wd)

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi_lite.Rds")

soi <- soi[, soi$seurat_clusters != "21"]
source("day9_only/step1.1a_helper_functions_2023-09-29.R")
soi <- t.f.cellType(
  soi, "day9_only/data/cell_type_correspondance_for_pct.tsv")

major.cellTypes <- c("HSPC", "Neu", "Mono/Mac/DC", "Ery", "Mega", "Lympho")

color.df <- read.table("day9_only/data/colorMap_for_pct.tsv",
                       comment.char = "", header = T, sep = "\t")
color.df <- color.df %>% filter(color != "")
colors <- color.df$color
names(colors) <- color.df$cell_type
colors <- colors[major.cellTypes]

soi$sample <- factor(soi$sample, levels = SAMPLES)

sample.map <- c("KO_d11_r1", "WT_d11_r1", "KO_d9_r1", "KO_d9_r2", "WT_d9_r1", "WT_d9_r2")
names(sample.map) <- c("KO_rep1", "WT_rep1", "KO_rep2", "KO_rep3", "WT_rep2", "WT_rep3")

soi$sample.pub <- sample.map[soi$sample]
soi$sample.pub <- factor(soi$sample.pub, levels = rev(sample.map))

p <- cluster.freq.bar(soi = soi, x = "sample.pub", group.by = "cell_type", groups = major.cellTypes,
                      average.across = "rep", average.method = "number", # not really useful since we are plotting each sample
                      color.map = colors,
                      out.file = paste0(wd, "/cluster_freq_majorClusters_bySample.pdf"))

#>>>>>>>>> spreadsheet for all clusters
soi <- t.f.cellType(
  soi, "day9_only/data/cellType_map.tsv")

count.by.cluster(so = soi, work.dir = wd, root.name = "cluster_freq", write.xls = T,
                 meta = "cell_type", group.by = "sample.pub", mixedSort = F)

