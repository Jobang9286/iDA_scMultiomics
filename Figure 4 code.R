#Advanced Science reproduce#

#Figure 4#

# read in peak sets
CONTROL <- read.table(
  file = "/home/student/13_scATAC-seq_IDA/CONTROL/BED/peaks.bed",
  col.names = c("chr", "start", "end")
)

D7 <- read.table(
  file = "/home/student/13_scATAC-seq_IDA/D7/BED/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
CONTROL <- makeGRangesFromDataFrame(CONTROL)
D7 <- makeGRangesFromDataFrame(D7)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(CONTROL, D7))

combined.peaks

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

md.con <- read.table(
  file = "/home/student/13_scATAC-seq_IDA/CONTROL/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.D7 <- read.table(
  file = "/home/student/13_scATAC-seq_IDA/D7/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.con <- md.con[md.con$passed_filters > 6000, ]
md.D7 <- md.D7[md.D7$passed_filters > 10000, ]

frags.con <- CreateFragmentObject(
  path = "/home/student/13_scATAC-seq_IDA/CONTROL/fragments.tsv.gz",
  cells = rownames(md.con)
)

frags.d7 <- CreateFragmentObject(
  path = "/home/student/13_scATAC-seq_IDA/D7/fragments.tsv.gz",
  cells = rownames(md.D7)
)

con.counts <- FeatureMatrix(
  fragments = frags.con,
  features = combined.peaks,
  cells = rownames(md.con)
)

d7.counts <- FeatureMatrix(
  fragments = frags.d7,
  features = combined.peaks,
  cells = rownames(md.D7)
)

con_assay <- CreateChromatinAssay(con.counts, fragments = frags.con)
con <- CreateSeuratObject(con_assay, assay = "ATAC", meta.data=md.con)

d7_assay <- CreateChromatinAssay(d7.counts, fragments = frags.d7)
D7 <- CreateSeuratObject(d7_assay, assay = "ATAC", meta.data=md.D7)

con$time <- 'MEF'
D7$time <- 'D7'

combined <- merge(
  x = con,
  y = list(D7),
  add.cell.ids = c("con", "D7")
)

combined$time <- factor(x = combined$time, levels = c("MEF", "D7"))

combined[["ATAC"]]
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:8)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:8)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)

combined_2 <- subset(combined, idents = c("Cluster_5"), invert = T)

DimPlot(combined_2, reduction = "umap", label = F, label.size = 5, pt.size = 0.8, cols = c("#46085C", "#3D4C8A", "#FDE725", "#D7E219")) #Figure 4A

DimPlot(combined_2, reduction = "umap", label = F, label.size = 5, pt.size = 0.8, group.by = "time", cols = c("#328491", "#E6B75D")) #Figure 4B

figure_4c_meta <- combined_2@meta.data
write.csv(figure_4c_meta, "./figure_4c_meta.csv") #Figure 4C

VlnPlot(combined_2, features = c("Acta2", "Foxa1", "Th", "Tagln", "Lmx1b", "Drd2", "Igfbp3", "Wnt5a", "Slc6a3"), pt.size = 0, ncol = 3, cols = c("#46085C", "#3D4C8A", "#FDE725", "#D7E219")) #Figure 4D

cluster_12_atac <- subset(combined_2, idents = c("Cluster_1", "Cluster_2"))

cluster_34_atac <- subset(combined_2, idents = c("Cluster_3", "Cluster_4"))

#Figure 4E

DefaultAssay(cluster_12_atac) <- 'ATAC'

active_meta <- cluster_12_atac@active.ident

cluster_12_atac <- AddMetaData(cluster_12_atac, metadata = active_meta, col.name = "seurat_clusters")

#obj_123_d <- subset(obj_123, cells = WhichCells(obj_123, downsample = 200)) down sampling하면 더 빠름

gene.activities <- GeneActivity(cluster_12_atac, features = VariableFeatures(obj_123_d))
cluster_12_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(cluster_12_atac) <- "ACTIVITY"
cluster_12_atac <- NormalizeData(cluster_12_atac)
cluster_12_atac <- ScaleData(cluster_12_atac, features = rownames(cluster_12_atac))

transfer.anchors <- FindTransferAnchors(reference = obj_123_d, query = cluster_12_atac, features = VariableFeatures(object = obj_123_d),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = obj_123_d$cell_type,
                                     weight.reduction = cluster_12_atac[["lsi"]], dims = 2:5)

cluster_12_atac <- AddMetaData(cluster_12_atac, metadata = celltype.predictions)

cluster_12_atac$annotation_correct <- cluster_12_atac$predicted.id == cluster_12_atac$seurat_clusters

DimPlot(cluster_12_atac, split.by = "predicted.id", group.by = "time") #Supple figure 3C

cluster_12_atac@meta.data

DimPlot(cluster_12_atac, pt.size = 0.9, label = T)

cluster_12_atac$seurat_clusters
predictions <- table(cluster_12_atac$seurat_clusters, cluster_12_atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var2, Var1, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#29788E") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p1 #Figure 4E



#Figure 4F

DefaultAssay(cluster_34_atac) <- 'ATAC'

#obj_45_d <- subset(obj_45, cells = WhichCells(obj_45, downsample = 200)) down sampling하면 더 빠름

gene.activities <- GeneActivity(cluster_34_atac, features = VariableFeatures(obj_45))
cluster_34_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(cluster_34_atac) <- "ACTIVITY"
cluster_34_atac <- NormalizeData(cluster_34_atac)
cluster_34_atac <- ScaleData(cluster_34_atac, features = rownames(cluster_34_atac))

transfer.anchors <- FindTransferAnchors(reference = obj_45, query = cluster_34_atac, features = VariableFeatures(object = obj_45),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = obj_45$cell_type,
                                     weight.reduction = cluster_34_atac[["lsi"]], dims = 2:5)


cluster_34_atac <- AddMetaData(cluster_34_atac, metadata = celltype.predictions)

cluster_34_atac$annotation_correct <- cluster_34_atac$predicted.id == cluster_34_atac$seurat_clusters

DimPlot(cluster_34_atac, split.by = "predicted.id", group.by = "time", pt.size = 0.7, cols = c("#29788E", "#E6B75B")) #Supple Figure 3D

cluster_34_atac$seurat_clusters

cluster_34_atac$seurat_clusters

predictions <- table(cluster_34_atac$seurat_clusters, cluster_34_atac$predicted.id)

predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type

predictions <- as.data.frame(predictions)

p2 <- ggplot(predictions, aes(Var2, Var1, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", low = "#ffffc8", high = "#29788E") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2 #Figure 4F

#Figure 4G
CoveragePlot(
  object = combined_2,
  region = c("Arid4b"),
  annotation = T,
  peaks = F
)

CoveragePlot(
  object = combined_2,
  region = c("Smarcb1"),
  annotation = T,
  peaks = F
)

CoveragePlot(
  object = combined_2,
  region = c("Smarcc2"),
  annotation = T,
  peaks = F
)

CoveragePlot(
  object = combined_2,
  region = c("Ascl1"),
  annotation = T,
  peaks = F
)

#Supple figure 4A

combined_2_d <- subset(combined_2, cells = WhichCells(combined_2, downsample = 200))

sce_1 <- as.SingleCellExperiment(combined_2_d, assay = "RNA")

shuffle <- sample(ncol(sce_1))

sce_1 <- slingshot(sce_1, reducedDim = 'UMAP', clusterLabels = sce_1$time, start.clus = 'Cluster_1', approx_points = 100)
layout(matrix(c(2, 2, 1, 1), 2))
par(mar = c(4.5, 4, 1, 1))

plot(reducedDims(sce_1)$UMAP, col = cell_colors, pch=16, asp = 1)

lines(SlingshotDataSet(sce_1), lwd=2, col='black')

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(sce_1$condition, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(sce_1$seurat_clusters, hue_pal())

plot(reducedDims(sce_1)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2",
     col = hcl.colors(100, alpha = .5)[cut(sce_1$slingPseudotime_1, breaks = 100)][shuffle])

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$time == "MEF", 1], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$time == "D7", 1], na.rm = T))


xlim <- range(c(ds$a_1$x, ds$a_2$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(10, "BrBG")[9], alpha = .8))
polygon(c(min(ds$a_2$x), ds$a_1$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(8, "Set2")[6], alpha = .8))

