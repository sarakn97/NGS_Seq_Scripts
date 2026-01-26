#### helper functions for seurat app

umap <- function(seu, grp, split, pt) {
  if (split == FALSE) {
    DimPlot(seu, group.by = grp, pt.size = pt)  +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))
  } else if (split == TRUE) {
    DimPlot(seu, group.by = grp, pt.size = pt, split.by = "sample") +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))
  }
}

featureplot <- function(seu, gene, lab, size, repell, point, color, split) {
  if (split == FALSE) {
    FeaturePlot(seu, feature = gene, label = lab, label.size = size, repel = repell, pt.size = point, cols = color, order = T) 
  } else if (split == TRUE) {
    FeaturePlot(seu, feature = gene, label = lab, label.size = size, repel = repell, pt.size = point, cols = color, order = T, split.by = "sample") 
  }
  
}

violin_plot <- function(seu, gene, split_vln) {
  if (split_vln == FALSE) {
    VlnPlot(seu, features = gene) + guides(fill = guide_legend(override.aes = list(size=4), ncol=1) )
  } else if (split_vln == TRUE) {
    VlnPlot(seu, features = gene,  split.by = "sample")
    }
}

# violin_plot <- function(seu, gene, group, split_vln) {
#   if (split_vln == FALSE) {
#     VlnPlot(seu, features = gene, group.by = group)
#   } else if (split_vln == TRUE) {
#     VlnPlot(seu, features = gene, group.by = group, split.by = "sample") }
# }