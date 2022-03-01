library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)

a<-c("G_N","G_T","G_T","G_N","G_T","G_N","G_T","G_T","G_N","G_T",
     "G_N","G_T","G_T","G_T","G_T","G_T","G_T","G_T","P_T","P_T",
     "G_N","G_T","G_N","G_T","G_N","G_T","G_T","G_T","G_T","G_T",
     "G_N","G_T","G_T","G_T","G_N","G_T","P_N","P_T","G_T","G_T")
for (i in 1:40) {
  assign(paste("sample_",i,"_data",sep = ""),read.csv(paste("D:/GSE183904_RAW/sample_",i,".csv",sep = ""),row.names = 1))
  b<-get(paste("sample_",i,"_data",sep = ""))
  for (j in 1:length(colnames(get(paste("sample_",i,"_data",sep = ""))))) {
    colnames(b)[j]<-paste(colnames(b[j]),a[i],j,sep = "-")
  }
  assign(paste("sample_",i,"_data",sep = ""),b)
  assign(paste("sample_",i,"_metadata",sep = ""),data.frame(colnames(get(paste("sample_",i,"_data",sep = ""))),rep(a[i],length(colnames(get(paste("sample_",i,"_data",sep = "")))))))
  c<-get(paste("sample_",i,"_metadata",sep = ""))
  d<-as.data.frame(colnames(get(paste("sample_",i,"_data",sep = ""))))
  colnames(c)<-c("barcode","group")
  rownames(c)<-d[,1]
}
e<-paste("sample_",1,"_metadata",sep="")

for (i in 2:40) {
  f<-paste("sample_",i,"_metadata",sep="") 
  e<-paste(e,f,sep = ", ")
}

sample_metadata<-rbind(sample_1_metadata, sample_2_metadata, sample_3_metadata, sample_4_metadata, sample_5_metadata, sample_6_metadata, sample_7_metadata, sample_8_metadata, sample_9_metadata, sample_10_metadata, sample_11_metadata, sample_12_metadata, sample_13_metadata, sample_14_metadata, sample_15_metadata, sample_16_metadata, sample_17_metadata, sample_18_metadata, sample_19_metadata, sample_20_metadata, sample_21_metadata, sample_22_metadata, sample_23_metadata, sample_24_metadata, sample_25_metadata, sample_26_metadata, sample_27_metadata, sample_28_metadata, sample_29_metadata, sample_30_metadata, sample_31_metadata, sample_32_metadata, sample_33_metadata, sample_34_metadata, sample_35_metadata, sample_36_metadata, sample_37_metadata, sample_38_metadata, sample_39_metadata, sample_40_metadata)
colnames(sample_metadata)<-c("barcodes","group")
rownames(sample_metadata)<-sample_metadata[,1]

remove(sample_1_metadata, sample_2_metadata, sample_3_metadata, sample_4_metadata, sample_5_metadata, sample_6_metadata, sample_7_metadata, sample_8_metadata, sample_9_metadata, sample_10_metadata, sample_11_metadata, sample_12_metadata, sample_13_metadata, sample_14_metadata, sample_15_metadata, sample_16_metadata, sample_17_metadata, sample_18_metadata, sample_19_metadata, sample_20_metadata, sample_21_metadata, sample_22_metadata, sample_23_metadata, sample_24_metadata, sample_25_metadata, sample_26_metadata, sample_27_metadata, sample_28_metadata, sample_29_metadata, sample_30_metadata, sample_31_metadata, sample_32_metadata, sample_33_metadata, sample_34_metadata, sample_35_metadata, sample_36_metadata, sample_37_metadata, sample_38_metadata, sample_39_metadata, sample_40_metadata,e,f)


e<-paste("sample_",1,"_data",sep="")

for (i in 2:40) {
  f<-paste("sample_",i,"_data",sep="") 
  e<-paste(e,f,sep = ", ")
}
sample_data<-cbind(sample_1_data, sample_2_data, sample_3_data, sample_4_data, sample_5_data, sample_6_data, sample_7_data, sample_8_data, sample_9_data, sample_10_data, sample_11_data, sample_12_data, sample_13_data, sample_14_data, sample_15_data, sample_16_data, sample_17_data, sample_18_data, sample_19_data, sample_20_data, sample_21_data, sample_22_data, sample_23_data, sample_24_data, sample_25_data, sample_26_data, sample_27_data, sample_28_data, sample_29_data, sample_30_data, sample_31_data, sample_32_data, sample_33_data, sample_34_data, sample_35_data, sample_36_data, sample_37_data, sample_38_data, sample_39_data, sample_40_data)

remove(sample_1_data, sample_2_data, sample_3_data, sample_4_data, sample_5_data, sample_6_data, sample_7_data, sample_8_data, sample_9_data, sample_10_data, sample_11_data, sample_12_data, sample_13_data, sample_14_data, sample_15_data, sample_16_data, sample_17_data, sample_18_data, sample_19_data, sample_20_data, sample_21_data, sample_22_data, sample_23_data, sample_24_data, sample_25_data, sample_26_data, sample_27_data, sample_28_data, sample_29_data, sample_30_data, sample_31_data, sample_32_data, sample_33_data, sample_34_data, sample_35_data, sample_36_data, sample_37_data, sample_38_data, sample_39_data, sample_40_data,e,f)
remove(a,b,c,d,i,j)

sample <- CreateSeuratObject(counts = sample_data, project = "sample",meta.data = sample_metadata,min.cells = 3, min.features = 200)
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by = "group")
sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 50)
sample <- NormalizeData(sample)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample)
sample <- ScaleData(sample, features = all.genes)
sample <- RunPCA(sample, features = VariableFeatures(object = sample))
ElbowPlot(sample)
sample <- FindNeighbors(sample, dims = 1:20)
sample <- FindClusters(sample, resolution = 1.2)
sample <- RunUMAP(sample, dims = 1:20)
sample <- RunTSNE(sample, dims = 1:20,check_duplicates = FALSE)
DimPlot(sample, reduction = "umap")
DimPlot(sample, reduction = "umap",split.by = "group")
DimPlot(sample, reduction = "tsne")
DimPlot(sample, reduction = "tsne",split.by = "group")
VlnPlot(sample,features = c("FOXP3"),pt.size = 0,sort = TRUE)
