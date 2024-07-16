#Advanced Science reproduce#

#Figure 2#

DimPlot(obj_1, reduction = "umap", pt.size = 0.7, label = F, label.size = 7, group.by = "condition", cols = c("#00BFC4", "#C578FF", "#7AAD00", "#F87269")) #figure 2A

obj_123 <- subset(obj_1, idents = c("Cluster_1","Cluster_2", "Cluster_3"))

obj_123 <- RunPCA(obj_123, npcs = 20, verbose = FALSE)
obj_123 <- RunUMAP(obj_123, reduction = "pca", dims = 1:5)
obj_123 <- FindNeighbors(obj_123, reduction = "pca", dims = 1:5)
obj_123 <- FindClusters(obj_123, resolution = 0.5)

DimPlot(obj_123, reduction = "umap", pt.size = 0.6, label = F, label.size = 7) #Figure 2B

hfib.genes <- as.matrix(read.csv("/home/student/11_DN_re/mouse scatter/hfib_as_mouse.csv"))
hfib.genes <- as.character(hfib.genes)

obj_123 = AddModuleScore(
  object = obj_123,
  features = list(hfib.genes),
  name = 'hfib.genes', ctrl = 3, assay = "RNA") #Femandes et al hFib

VlnPlot(obj_123, features = "hfib.genes1", pt.size = 0) +
  geom_boxplot(outlier.size = 0) #Figure 2C

#DEG참고내용 : adjust pvalue 0은 제외한 logFC 0.25 이상 / -0.25 이하 유전자로 DEG number를 산정하였음

C1_vs_C2 <- FindMarkers(obj_123, ident.1 = "Cluster_1", ident.2 = "Cluster_2", logfc.threshold = 0.1) #Figure 2D and 2E
C2_vs_C3 <- FindMarkers(obj_123, ident.1 = "Cluster_2", ident.2 = "Cluster_3", logfc.threshold = 0.1) #Figure 2D
C3_vs_C1 <- FindMarkers(obj_123, ident.1 = "Cluster_3", ident.2 = "Cluster_1", logfc.threshold = 0.1) #Fiutre 2D and 2F


#Figure 2G

sce_1 <- as.SingleCellExperiment(obj_123, assay = "RNA")

shuffle <- sample(ncol(sce_1))

sce_1 <- slingshot(sce_1, reducedDim = 'UMAP', clusterLabels = sce_1$ident, start.clus = 'Cluster_1', approx_points = 100)
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

ds <- list(a_1 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "fib", 1], na.rm = T),
           a_2 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d7_obj", 1], na.rm = T),
           a_3 = density(slingPseudotime(sce_1)[colData(sce_1)$condition == "dn_d10_obj", 1], na.rm = T))

xlim <- range(c(ds$a_1$x, ds$a_2$x))
ylim <- range(c(ds$a_1$y, ds$a_2$y))

plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")

polygon(c(min(ds$a_1$x), ds$a_1$x, max(ds$a_1$x)), c(0, ds$a_1$y, 0),
        col = alpha(brewer.pal(10, "Set3")[9], alpha = .5))
polygon(c(min(ds$a_2$x), ds$a_2$x, max(ds$a_2$x)), c(0, ds$a_2$y, 0),
        col = alpha(brewer.pal(10, "Set3")[1], alpha = .5))
polygon(c(min(ds$a_3$x), ds$a_3$x, max(ds$a_3$x)), c(0, ds$a_3$y, 0),
        col = alpha(brewer.pal(10, "Reds")[8], alpha = .5)) #색 코드 까먹음

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

plotSmoothers(sce, assays(sce)$counts, gene = c("Hoxa10"), alpha = 0.6, border = T, lwd = 2) #Figure 2J
plotSmoothers(sce, assays(sce)$counts, gene = c("Hoxa11os"), alpha = 0.6, border = T, lwd = 2)
plotSmoothers(sce, assays(sce)$counts, gene = c("Hoxa13"), alpha = 0.6, border = T, lwd = 2)
plotSmoothers(sce, assays(sce)$counts, gene = c("Hoxd13"), alpha = 0.6, border = T, lwd = 2)



