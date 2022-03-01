setwd("E:/AM20025-1/YL_analysis/homer_WT_KO")
a = list.files("peaks_regions")   
dir = paste("./peaks_regions/",a,sep="")
n = length(dir)  
c = read.table("E:/AM20025-1/YL_analysis/homer_WT_KO/column.txt",sep = "\t")
for (i in seq(1,n,2) ) {
  merge.data = read.table(file = dir[i],sep="\t")
  merge.data2 = read.table(file = dir[i+1],sep="\t",header = TRUE)
  merge.data_order <- merge.data[order(merge.data[,1]), ]
  merge.data_order2 <- merge.data2[order(merge.data2[,1]), ]
  merge.new <-cbind(merge.data_order,merge.data_order2)
  colnames(merge.new)<-c[1, ]
  write.table(merge.new,file =a[i],sep="\t",row.names = FALSE)
}
