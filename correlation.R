library(dplyr)
library(Seurat)
options(warn=-1)
pbmc.data <- Read10X(data.dir="I:/BOAO single-cell sequencing/BOAO/GCB1/2.2.filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB1_filtered", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:15, verbose=F)
pbmc <- FindClusters(pbmc, resolution = 0.6, verbose=F)
VlnPlot(pbmc,features = c("Kdm6b","Myc","Irf4","Prdm1"))
pbmc <- RunTSNE(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "tsne")

GCB1_nomrlize_scaledata<-pbmc@assays[["RNA"]]@scale.data
GCB1_nomrlize_scaledata<-as.data.frame(GCB1_nomrlize_scaledata)

GCB1_nomrlize_scaledata_t<-t(GCB1_nomrlize_scaledata)
GCB1_nomrlize_scaledata_t<-as.data.frame(GCB1_nomrlize_scaledata_t)

kdm6b_Nur77<-filter(GCB1_nomrlize_scaledata_t,Kdm6b >0 & Nr4a1 >0)
R<-rownames(kdm6b_Nur77)

library(ggplot2)
library(ggpubr)
ggplot(data = kdm6b_Nur77,aes(x=Kdm6b,y=Nr4a1)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = kdm6b_Nur77,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
FeaturePlot(pbmc,features = c("Kdm6b","Nr4a1"),cells = R,label = TRUE,pt.size = 3.4)


