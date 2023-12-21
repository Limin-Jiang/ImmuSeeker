library(data.table)
library(data.table)
library(stringr)
library(docopt)

doc <- "
Usage: myscript.R [-i <Imputname>] [-e <error>] [-c <threshold>] [-p <probability>] [-x <expression>]

Options:
  -i, --Imputname   The file of list of cancers
  -e, --error       Specify the sequencing error ratio (default: 0.02)
  -c, --threshold   Specify min number of supported reads (default: 0)
  -p, --probability Specify the probability of three genotypes with '-p (1/4,1/2,1/4)' (default: (1/3,1/3,1/3))
  -x, --expression  Specify whether to detect the HLA gene expression value
"

options <- docopt(doc)

## parameters
file_ID = options$Imputname #file_ID = "output.Nuc.txt_All.csv"
erp = as.numeric(options$error) #erp = 0.02
erp = 1 - erp
threshold = as.numeric(options$threshold) #threshold = 0
Pp = options$probability #Pp = "(1/3,1/3,1/3)"
Pp = eval(parse(text = paste0("c", Pp)))
exlag = options$expression
#print(erp)
#print(threshold)
#print(Pp)

#setwd("C:\\Users\\lxj423\\OneDrive - University of Miami\\work_data\\HLA\\NUC_nEW\\R\\pipeline")

getP = function(data_temp,Pp,erp){
  numsum = sum(data_temp$Num )
  binomialcoff = data_temp$Num[1]
  Pd = c(dbinom(binomialcoff, size=numsum, prob=erp), dbinom(binomialcoff, size=numsum, prob=0.5),dbinom(binomialcoff, size=numsum, prob=(1-erp) )  )
  
  PAAd = (Pd[1]* Pp[1])/(sum(Pd * Pp))
  PAGd = (Pd[2]* Pp[2])/(sum(Pd * Pp))
  PGGd = (Pd[3]* Pp[3])/(sum(Pd * Pp))
  
  data_end = data.frame(Genetype = c(paste0(data_temp$HLAname[1],",",data_temp$HLAname[1]),paste0(data_temp$HLAname[1],",",data_temp$HLAname[2]),paste0(data_temp$HLAname[2],",",data_temp$HLAname[2])),Pro = c(PAAd,PAGd,PGGd),
                        Number = c(paste0(data_temp$Num[1],"?",data_temp$Num[1]),paste0(data_temp$Num[1],"?",data_temp$Num[2]),paste0(data_temp$Num[2],"?",data_temp$Num[2])),stringsAsFactors = FALSE)
  data_end$PL = -10*log10(data_end$Pro + 10^-100)
  data_end$Normalized_PL = data_end$PL - min(data_end$PL)
  data_end$Type = c("homo","heter","homo")
  return(data_end)
}


Getdata = function(subdata){
  
  if (max(subdata$Num) > 1000){
    subdata$Num =round(subdata$Num * (1000/max(subdata$Num)))
  }
  
  
  if (nrow(subdata) > 1){
    subdata = subdata[order(subdata$Num,decreasing=TRUE),]
    data_temp = subdata[c(1:2),]
    data_temp = getP(data_temp,Pp,erp)
    data_temp$gene = gene_ID
    data_temp$Sample = basename(file_ID) 
  }
  if (nrow(subdata) == 1){
    subdata = subdata[order(subdata$Num,decreasing=TRUE),]
    data_temp = rbind(subdata,subdata)
    data_temp$HLAname[2]= ""
    data_temp$Num[2] = 0
    data_temp = getP(data_temp,Pp,erp)
    data_temp$gene = gene_ID
    data_temp$Sample = basename(file_ID) 
    data_temp = data_temp[1,]
  }
  return(data_temp)
  
}



data_all = data.frame()
d1 = fread(file_ID,header = T,sep = ",")
d1 = subset(d1,Num>threshold)
geness = unique(d1$Genes)
data_temp = data.frame()

for (gene_ID in geness){
  subdata = subset(d1, Genes == gene_ID & Class == "Class1" & Num > 0)
  if (nrow(subdata) > 0){
    data_all = rbind(data_all,Getdata(subdata))
  }
}


data_all$Sample = gsub("_HLAlist.csv","",data_all$Sample)



d2 = data_all

d2$ID = paste(d2$Sample,d2$gene,sep = "_")
IDs = unique(d2$ID)


data_all = data.frame()
for (IDs_id in IDs){
  subd1 = subset(d2,ID == IDs_id)
  subd1 = subd1[which.max(subd1$Pro),]
  if (nrow(subd1) > 0){
    data_all = rbind(data_all,subd1)
  }
}


data_all = data_all[,-ncol(data_all)]
write.csv(data_all, paste0(gsub("_HLAlist.csv","",file_ID),"_genetype.csv") ,row.names = FALSE)

if (exlag == "true"){
	Hla_lists = unique(as.vector(sapply(strsplit(data_all$Genetype, split=',', fixed=TRUE), function(x) (x[]))))
	subd1 = subset(d1,Class == "Class1" & HLAname %in% Hla_lists)
	d_search = subd1
	subd1 = subset(d1,Class == "Class2" & Type == "Unique")
	subd1$lags = sapply(strsplit(subd1$HLAname, split=':', fixed=TRUE), function(x) (x[1]))
	subd1 = subset(subd1,lags %in% Hla_lists)
	subd1 = subd1[,-6]
	d_search = rbind(d_search,subd1)

	write.csv(d_search, paste0(file_ID,".search") ,row.names = FALSE)
	d_search = unique(d_search$Genes)
	d_search = data.frame(Gene_ID = d_search, Gene_ID1 = d_search, stringsAsFactors = FALSE )
	write.csv(d_search, paste0(file_ID,".gene.search") ,row.names = FALSE)
}

