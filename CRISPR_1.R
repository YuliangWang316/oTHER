setwd("E:/GSE127232/")
a = list.files("RAW")   
n = length(a)
e = "mageck count -l library.csv -n result --sample-label"
for (i in 5:n) {
  c = strsplit(a[i],'_')
  e = paste(e,c[[1]][3],sep = ",")
}
d = paste(e," --fastq")
for (k in 5:n) {
  d = paste(d,a[k],sep = " ")
}
d
