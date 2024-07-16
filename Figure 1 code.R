#Advanced Science reproduce#

#Figure 1#

library(SAVER)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(cowplot)
library(slingshot)
library(tradeSeq)
library(slingshot)
library(Scillus)
library(tidyverse)
library(magrittr)
library(clusterProfiler)
library(gridExtra)

setwd("/home/student/31_iDA_reproduce/")

fib <- Seurat::Read10X("/home/student/11_DN_re/WT/outs/filtered_feature_bc_matrix/")
fib_obj <- CreateSeuratObject(fib, project = "fib", min.cells = 3, min.features = 200)
fib_obj$condition <- "fib"
fib_obj$sample <- "JY"
fib_obj[["percent.mt"]] <- PercentageFeatureSet(fib_obj, pattern = "^mt-")
VlnPlot(fib_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fib_obj <- subset(fib_obj, subset = nFeature_RNA > 200 & nCount_RNA < 40000 & nCount_RNA > 1000 & percent.mt < 10)

dn_d7 <- Seurat::Read10X("/home/student/11_DN_re/D7/outs/filtered_feature_bc_matrix/")
dn_d7_obj <- CreateSeuratObject(dn_d5, project = "dn_d7", min.cells = 3, min.features = 200)
dn_d7_obj$condition <- "dn_d7_obj"
dn_d7_obj$sample <- "JY"
dn_d7_obj[["percent.mt"]] <- PercentageFeatureSet(dn_d7_obj, pattern = "^mt-")
VlnPlot(dn_d7_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dn_d7_obj <- subset(dn_d7_obj, subset = nFeature_RNA > 200 & nCount_RNA < 70000 & nCount_RNA > 2000 & percent.mt < 10)

dn_d10 <- Seurat::Read10X("/home/student/11_DN_re/D10/outs/filtered_feature_bc_matrix/")
dn_d10_obj <- CreateSeuratObject(dn_d11,project = "dn_d10", min.cells = 3, min.features = 200)
dn_d10_obj$condition <- "dn_d10_obj"
dn_d10_obj$sample <- "JY"
dn_d10_obj[["percent.mt"]] <- PercentageFeatureSet(dn_d10_obj, pattern = "^mt-")
VlnPlot(dn_d10_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dn_d10_obj <- subset(dn_d10_obj, subset = nFeature_RNA > 200 & nCount_RNA < 30000 & nCount_RNA > 1000 & percent.mt < 10)

dn_d21 <- Seurat::Read10X("/home/student/11_DN_re/D21/outs/filtered_feature_bc_matrix/")
dn_d21_obj <- CreateSeuratObject(dn_d30,project = "dn_d21", min.cells = 3, min.features = 200)
dn_d21_obj$condition <- "dn_d21_obj"
dn_d21_obj$sample <- "JY"
dn_d21_obj[["percent.mt"]] <- PercentageFeatureSet(dn_d21_obj, pattern = "^mt-")
VlnPlot(dn_d21_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dn_d21_obj <- subset(dn_d21_obj, subset = nFeature_RNA > 200 & nCount_RNA < 50000 & nCount_RNA > 3000 & percent.mt < 10)


obj <- merge(x= fib_obj, y=list(dn_d7_obj, dn_d10_obj, dn_d21_obj))
view(obj@meta.data)
class(obj[["RNA"]])

obj <- JoinLayers(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj
class(obj[["RNA"]])
Layers(obj[["RNA"]])

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj)
obj <- FindClusters(obj, resolution = 0.8)
obj <- RunUMAP(obj, dims = 1:12)
DimPlot(obj)
DimPlot(obj, reduction = "umap", pt.size = 0.7, label = F, label.size = 7, cols = c("#480C5E", "#3D4D8A", "#328293", "#77D574", "#CCE12C", "#FDE725")) #figure 1B upper

meta_obj <- obj@meta.data
meta_obj_1 <- obj@active.ident
meta_obj <- cbind(meta_obj, meta_obj_1)
meta_obj
write.csv(meta_obj, "./meta_obj.csv") #figure 1B bottom

obj_down <- subset(obj, cells = WhichCells(obj, downsample = 300))


obj.markers <- FindAllMarkers(obj_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC) -> top10

top10$gene

DoHeatmap(obj_down, features = top10$gene)  + scale_fill_gradientn(colors = c("#46085C", "#23A488", "#FDE622")) #figure 1C

hfib.genes <- as.matrix(read.csv("/home/student/11_DN_re/mouse scatter/hfib_as_mouse.csv"))
iN.genes <- as.matrix(read.csv("/home/student/11_DN_re/iN.csv"))

hDA1.genes <- as.matrix(read.csv("/home/student/11_DN_re/mouse scatter/hDA_1_as.csv"))
hDA4.genes <- as.matrix(read.csv("/home/student/11_DN_re/mouse scatter/hDA_4_as.csv"))

hfib.genes <- as.character(hfib.genes)
iN.genes <- as.character(iN.genes)

hDA1.genes <- as.character(hDA1.genes)
hDA4.genes <- as.character(hDA4.genes)

obj = AddModuleScore(
  object = obj,
  features = list(hfib.genes),
  name = 'hfib.genes', ctrl = 3, assay = "RNA") #Femandes et al hFib

obj = AddModuleScore(
  object = obj,
  features = list(iN.genes),
  name = 'iN.genes', ctrl = 3, assay = "RNA") #Femandes et al iN

obj = AddModuleScore(
  object = obj,
  features = list(hDA1.genes),
  name = 'hDA1.genes', ctrl = 3, assay = "RNA") #Cates et al DN 1

obj = AddModuleScore(
  object = obj,
  features = list(hDA4.genes),
  name = 'hDA4.genes', ctrl = 3, assay = "RNA") #Cates et al DN 2

features <- c("Tagln", "Acta2", "Ogn", "Mki67", "Ccnb2", "Pbx1", "Gap43", "Stmn2",
              "Mnat", "Gng3", "Srxn1", "Gclc")

VlnPlot(obj, features, stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("Representative markers of Cell type") #Figure 1D

DotPlot(obj, features = c("hfib.genes1", "iN.genes1","hDA1.genes1", "hDA4.genes1"), cols = c("yellow", "red")) #figure 1E

FeatureScatter(obj, feature1 = "hDA1.genes1", feature2 = "hfib.genes1", group.by = "condition", cols = c("#00BFC4", "#C578FF", "#7AAD00", "#F87269")) #Figure 1F





                                                              