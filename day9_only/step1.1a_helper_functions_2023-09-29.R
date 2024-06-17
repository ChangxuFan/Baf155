t.f.cellType <- function(so, df, from = "seurat_clusters") {
  # df: 2 columns, corresponding to source and target clusters. 
  # Thus forming the correspondance between them
  if (is.character(df)) {
    df <- read.table(df, header = T)
    df[, from] <- as.character(df[, from])
  }
  
  if (any(duplicated(df[, from]))) {
    stop("some of the clusters in df$from is duplicated")
  }
  
  to <- colnames(df) %>% .[. != from]
  
  meta <- so@meta.data
  meta$cells <- rownames(meta)
  
  meta[[to]] <- NULL
  
  if (is.factor(meta[, from])) {
    levels <- levels(meta[, from])
    meta[, from] <- as.character(meta[, from])
  }
  meta2 <- dplyr::left_join(meta, df, by = from)
  
  if (!identical(meta$cells, meta2$cells)) {
    stop("something is wrong")
  }
  
  rownames(meta2) <- meta2$cells
  
  meta2[, from] <- factor(meta2[, from], levels = levels)
  
  so@meta.data <- meta2
  return(so)
}


t.f.plot.umap <- function(soi, aoi, rna = T, atac = T, out.dir, cell.type.map, color.map) {
  if (!rna && !atac) {
    stop("nothing to do")
  }
  if (any(soi$seurat_clusters == "21")) {
    stop("seurat cluster 21 should have been removed from soi")
  }
  
  aoi <- archr.add.seurat(ao = aoi, so = soi, meta = "cell_type", as.factor = F)
  aoi <- archr.add.seurat(ao = aoi, so = soi, meta = "sample", as.factor = F)
  aoi <- aoi[aoi$seurat_clusters %in% as.character(0:13),]
  
  if (any(aoi$seurat_clusters %in% as.character(14:21)))   {
    stop("small clusters (starting from 14) should have been removed from aoi")
  }
  
  soi <- t.f.cellType(soi, cell.type.map)
  
  color.df <- read.table(color.map, comment.char = "", header = T, sep = "\t")
  color.df <- color.df %>% filter(color != "")
  
  colors <- color.df$color
  names(colors) <- color.df$cell_type
  
  lapply(c("rna", "atac"), function(dataType) {
    if (dataType == "atac" && !atac)
      return()
    if (dataType == "rna" && !rna)
      return()
    
    lapply(c("d9", "d11", "all"), function(day) {
      samples <- "all"
      if (day == "d9") {
        samples <- c("KO_rep2", "KO_rep3", "WT_rep2", "WT_rep3")
      } else if (day == "d11") {
        samples <- c(c("KO_rep1", "WT_rep1"))
      } 
      
      if (samples != "all") {
        soi <- soi[, soi$sample %in% samples]
        aoi <- aoi[aoi$Sample %in% samples, ]
      }
      
      lapply(c("cell_type", "sample"), function(ident) {
        out.file <- paste0(out.dir, "/", dataType, "_", day, "_", ident, ".pdf")
        print(paste0(">>> preparing plot for: ", out.file))
        
        if (ident == "sample") {
          colors <- NULL
        }
        if (dataType == "rna") {
          umap.fanc(obj = soi, group.by = ident, cols = colors, pt.size = 0.2, 
                    plot.out = out.file)
        } else {
          umap.fanc(obj = aoi, group.by = ident, cols = colors, pt.size = 0.35, 
                    plot.out = out.file)
        }
        
      })
    })
  })
  
}
