#Advanced Science reproduce#

#Figure 3#

obj_456 <- subset(obj_1, idents = c("Cluster_4","Cluster_5", "Cluster_6"))
obj_456 <- RunPCA(obj_456, npcs = 20, verbose = FALSE)
obj_456 <- RunUMAP(obj_456, reduction = "pca", dims = 1:15)
obj_456 <- FindNeighbors(obj_456, reduction = "pca", dims = 1:15)
obj_456 <- FindClusters(obj_456, resolution = 0.5)

DimPlot(obj_456, reduction = "umap", pt.size = 0.6, label = F, label.size = 7) #Figure 3A

#Figure 3B-C

sce_1 <- as.SingleCellExperiment(obj_456, assay = "RNA")

shuffle <- sample(ncol(sce_1))

sce_1 <- slingshot(sce_1, reducedDim = 'UMAP', clusterLabels = sce_1$ident, start.clus = 'Cluster_4', approx_points = 100)
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
     col = hcl.colors(100, alpha = .5)[cut(sce_1$slingPseudotime_1, breaks = 100)][shuffle]) #Figure 3C

plot(reducedDims(sce_1)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2",
     col = hcl.colors(100, alpha = .5)[cut(sce_1$slingPseudotime_2, breaks = 100)][shuffle]) #Figure 3B

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d7_obj", 2], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d10_obj", 2], na.rm = T),
           a_3 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d21_obj", 2], na.rm = T)) #Figure 3B

xlim <- range(c(ds$a_1$x, ds$a_2$x, ds$a_3$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y, ds$a_3$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(10, "Set3")[9], alpha = .5))
polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(10, "Set3")[1], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(10, "Reds")[8], alpha = .5)) #Figure 3B #색코드 까먹음

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d7_obj", 1], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d10_obj", 1], na.rm = T),
           a_3 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d21_obj", 1], na.rm = T)) #Figure 3C

xlim <- range(c(ds$a_1$x, ds$a_2$x, ds$a_3$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y, ds$a_3$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(10, "Set3")[9], alpha = .5))
polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(10, "Set3")[1], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(10, "Reds")[8], alpha = .5)) #Figure 3C #색코드 까먹음

#Figure 3D

sce <- sce_1

slingsce<-SlingshotDataSet(sce)

pseudotimeED <- slingPseudotime(slingsce, na = FALSE)

cellWeightsED <- slingCurveWeights(slingsce)
conditions <- colData(sce)$condition
counts<-sce@assays@data@listData$counts
set.seed(3)

sce <- fitGAM(counts = counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T, genes = 1:2500, conditions = as.factor(sce$condition))

mean(rowData(sce)$tradeSeq$converged)

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))

assocRes <- rowData(sce)$assocRes
assocRes

plotSmoothers(sce, assays(sce)$counts, gene = c("Otx2"), alpha = 0.6, border = T, lwd = 2) #<- Otx2만 fitgam gene 수 늘려야 나오고 밑에 8개 gene들은 도출가능
plotSmoothers(sce, assays(sce)$counts, gene = c("En2"), alpha = 0.6, border = T, lwd = 2)#Figure 3D
plotSmoothers(sce, assays(sce)$counts, gene = c("Bcl11a"), alpha = 0.6, border = T, lwd = 2)#Figure 3D
plotSmoothers(sce, assays(sce)$counts, gene = c("Col3a1"), alpha = 0.6, border = T, lwd = 2)#Figure 3D
plotSmoothers(sce, assays(sce)$counts, gene = c("Col1a1"), alpha = 0.6, border = T, lwd = 2)#Figure 3D
plotSmoothers(sce, assays(sce)$counts, gene = c("Igfbp5"), alpha = 0.6, border = T, lwd = 2)#Figure 3D
plotSmoothers(sce, assays(sce)$counts, gene = c("Ptgs2"), alpha = 0.6, border = T, lwd = 2)#Figure 3D
plotSmoothers(sce, assays(sce)$counts, gene = c("S100a8"), alpha = 0.6, border = T, lwd = 2)#Figure 3D
plotSmoothers(sce, assays(sce)$counts, gene = c("Perp"), alpha = 0.6, border = T, lwd = 2)#Figure 3D

C5_vs_C6 <- FindMarkers(obj_456, ident.1 = "Cluster5", ident.2 = "Cluster_6", logfc.threshold = 0.1) #Figure 3H

VlnPlot(obj_456, features = c("Arid4b", "Srxn1", "Smarcc2", "Hmox1", "Smarcb1", "Areg"), pt.size = 0, ncol = 2) #Figure 3I

plotSmoothers(sce, assays(sce)$counts, gene = c("Arid4b"), alpha = 0.6, border = T, lwd = 2)#Figure 3J
plotSmoothers(sce, assays(sce)$counts, gene = c("Smarcc2"), alpha = 0.6, border = T, lwd = 2)#Figure 3J
plotSmoothers(sce, assays(sce)$counts, gene = c("Smarcb1"), alpha = 0.6, border = T, lwd = 2)#Figure 3J
plotSmoothers(sce, assays(sce)$counts, gene = c("Srxn1"), alpha = 0.6, border = T, lwd = 2)#Figure 3J
plotSmoothers(sce, assays(sce)$counts, gene = c("Hmox1"), alpha = 0.6, border = T, lwd = 2)#Figure 3J
plotSmoothers(sce, assays(sce)$counts, gene = c("Areg"), alpha = 0.6, border = T, lwd = 2)#Figure 3J
