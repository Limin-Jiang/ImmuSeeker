library(data.table)
library(stringr)
library(docopt)

doc <- "
Usage: myscript.R [-i <Imputname>]

Options:
  -i, --Imputname   The file of list of cancers
"
options <- docopt(doc)

## parameters
file_ID = options$Imputname #file_ID = "output"

#setwd("C:\\Users\\lxj423\\OneDrive - University of Miami\\work_data\\HLA\\NUC_nEW\\R\\pipeline")

d_total = fread(paste0(file_ID,".Reads.number.total"),header = FALSE,sep = " ")
d_total = d_total$V2/10^6

d_gene_length = fread("data/HLA-length.txt",header = FALSE,sep = "\t")
colnames(d_gene_length) = c("HLAID","HLAname","length")
d_gene_length$length = d_gene_length$length/1000


d1_Count = fread(paste0(file_ID,"_HLAlist.csv.search1.count"),header = FALSE,sep = " ")
d1_Count1 = fread(paste0(file_ID,"_HLAlist.csv.gene.search1.count"),header = FALSE,sep = " ")
d1_Count = rbind(d1_Count,d1_Count1)
d1_Count$V1 = gsub("'","",d1_Count$V1)
d1_Count = subset(d1_Count,V2>0)
colnames(d1_Count)[1:2] = c("HLAname","ReadsNum")


d1_Count = merge(d1_Count,d_gene_length,by.x = "HLAname",by.y = "HLAname",all.x = TRUE)

for (ii in 1:nrow(d1_Count)){
  HLA_idx = d1_Count$HLAname[ii]
  d_gene_length_temp = d_gene_length[grepl(paste0("^",gsub("\\*","\\\\*",HLA_idx)),d_gene_length$HLAname) ,]  
  d1_Count$length[ii] = median(d_gene_length_temp$length)
}


d1_Count$RPK = d1_Count$ReadsNum/d1_Count$length
d1_Count$RPKM = (d1_Count$RPK*10^6)/rep(d_total,nrow(d1_Count))
d1_Count$TPM = (d1_Count$RPKM*10^6)/(sum(d1_Count$RPKM,na.rm = TRUE))


d1_Count$totalNum = rep(d_total*10^6,nrow(d1_Count))

d1_Count$Level = rep("Gene",nrow(d1_Count))
ll_idx = which(str_count(d1_Count$HLAname, ":") == 1)
if (length(ll_idx) > 0){
  d1_Count$Level[ll_idx] = rep("2",length(ll_idx))
}
ll_idx = which(str_count(d1_Count$HLAname, ":") == 0 & grepl("\\*",d1_Count$HLAname))

if (length(ll_idx) > 0){
  d1_Count$Level[ll_idx] = rep("1",length(ll_idx))
}
d1_Count$length = d1_Count$length * 1000


idx = which(d1_Count$Level == "1")

subd1_Count = d1_Count[idx,]
d1_Count = d1_Count[-idx,]
d1_Count$Type = rep("NA",nrow(d1_Count))


d_homo = fread(paste0(file_ID,"_genetype.csv"),header = TRUE,sep = ",")
d_homo = d_homo[,c("gene", "Type")]
subd1_Count$gene = sapply(strsplit(subd1_Count$HLAname, split='*', fixed=TRUE), function(x) (x[1]))
subd1_Count = merge(subd1_Count,d_homo,by = "gene")
subd1_Count = subd1_Count[,-1]

d1_Count = rbind(subd1_Count,d1_Count)

d1_Count = d1_Count[order(d1_Count$Level),]


write.csv(d1_Count,paste0(file_ID,".expression.csv"),row.names = FALSE)

