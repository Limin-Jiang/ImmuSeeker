rm(list = ls())
library(data.table)
library(stringr)
#library(dplyr)
args=commandArgs(TRUE)
file_ID = args[1] 
print(file_ID)

#file_ID = "test1.txt"
#setwd("C:\\Users\\lxj423\\OneDrive - University of Miami\\work_data\\HLA\\Alignment\\match_all")


Newbam = fread(file_ID,header = F,sep = " ")
colnames(Newbam)[1] = "reads"

geneNameS = unique(sapply(strsplit(unique(Newbam$V2), split='*', fixed=TRUE), function(x) (x[1])))



data_all = data.frame()

for (geneName in geneNameS){
  
  
  print(geneName)
  
  subd1_temp = Newbam[!grepl(paste0("^",geneName,"\\*"), Newbam$V2),]
  subd1_temp = unique(subd1_temp$reads)
  subd1 = Newbam[grepl(paste0("^",geneName,"\\*"), Newbam$V2),]
  subd1 = subset(subd1,!(reads %in% subd1_temp))
  
  if (nrow(subd1) > 0){
    
    aa = sapply(strsplit(subd1$V2, split=':', fixed=TRUE), function(x) (x[1:4]))
    
    subd1$Class1 = aa[1,]
    subd1$Class2 = paste(aa[1,],aa[2,],sep = ":")
    subd1$Class3 = paste(aa[1,],aa[2,],aa[3,],sep = ":")
    
    colnames(subd1)[2] = "Class4"
    
    
    subd1$Class1[grepl(":NA",subd1$Class1)] = NA
    subd1$Class2[grepl(":NA",subd1$Class2)] = NA
    subd1$Class3[grepl(":NA",subd1$Class3)] = NA
    subd1$Class4[grepl(":NA",subd1$Class4)] = NA
    
    data_temp4 = data.frame(table(subd1$Class4))
    data_temp4$Var1 =  as.character(data_temp4$Var1)
    data_temp3 = data.frame(table(subd1$Class3))
    data_temp3$Var1 =  as.character(data_temp3$Var1)
    data_temp2 = data.frame(table(subd1$Class2))
    data_temp2$Var1 =  as.character(data_temp2$Var1)
    data_temp1 = data.frame(table(subd1$Class1))
    data_temp1$Var1 =  as.character(data_temp1$Var1)
    
    Class4_list = data_temp4$Var1[str_count(data_temp4$Var1, ":") == 3]
    Class3_list = data_temp3$Var1[str_count(data_temp3$Var1, ":") == 2]
    Class2_list = data_temp2$Var1[str_count(data_temp2$Var1, ":") == 1]
    Class1_list = data_temp1$Var1[str_count(data_temp1$Var1, ":") == 0]
    
    
    data_temp = data.frame(HLAname = geneName, Num = length(unique(subd1$reads)),Class = "Gene", Genes = geneName, Type = "Number", stringsAsFactors = FALSE)
    data_all = rbind(data_all,data_temp)
    
    
    for (Class1_list_ID in Class1_list){
      subsubset1 = subset(subd1, Class1 != Class1_list_ID)
      subsubset1 = unique(subsubset1$reads)
      subsubd1 = subset(subd1, Class1 == Class1_list_ID & !(reads %in% subsubset1))
      
      Class4_list_temp = Class4_list[grepl(paste0("^",  gsub("\\*","\\\\*",Class1_list_ID),":" ),Class4_list)]
      Class3_list_temp = Class3_list[grepl(paste0("^",  gsub("\\*","\\\\*",Class1_list_ID),":" ),Class3_list)]
      Class2_list_temp = Class2_list[grepl(paste0("^",  gsub("\\*","\\\\*",Class1_list_ID),":" ),Class2_list)]
      
      data_temp = data.frame(HLAname = Class1_list_ID, Num =  length(unique(subsubd1$reads )), Class = "Class1", Genes = geneName, Type = "Unique", stringsAsFactors = FALSE)
      data_all = rbind(data_all,data_temp)
      
      
      if (nrow(subsubd1) > 0){
        
        data_temp = data.frame()
        for (ii in Class2_list_temp){
          subsubset2 = subset(subsubd1, Class2 == ii )
          aa = unique(subsubset2$reads )
          data_temp = rbind(data_temp,data.frame(HLAname = ii, Num = length(aa),stringsAsFactors = FALSE))
        }
        
        if (nrow(data_temp) > 0){
          data_temp$Class = "Class2"
          data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
          data_temp$Genes = geneName
          data_temp$Type = "Number"
          data_all = rbind(data_all,data_temp)
        }
        
        
        
        
        data_temp = data.frame()
        for (ii in Class2_list_temp){
          subsubset1 = subset(subsubd1, Class2 != ii )
          subsubset2 = subset(subsubd1, Class2 == ii )
          aa = setdiff(subsubset2$reads,subsubset1$reads)
          data_temp = rbind(data_temp,data.frame(HLAname = ii, Num = length(aa),stringsAsFactors = FALSE))
          
        }
        if (nrow(data_temp) > 0){
          data_temp$Class = "Class2"
          data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
          data_temp$Genes = geneName
          data_temp$Type = "Unique"
          data_all = rbind(data_all,data_temp)
        }
        
        
        
        
        
        data_temp = data.frame()
        for (ii in Class3_list_temp){
          subsubset2 = subset(subsubd1, Class3 == ii )
          aa = unique(subsubset2$reads )
          data_temp = rbind(data_temp,data.frame(HLAname = ii, Num = length(aa),stringsAsFactors = FALSE))
          
        }
        
        if (nrow(data_temp) > 0){
          data_temp$Class = "Class3"
          data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
          data_temp$Genes = geneName
          data_temp$Type = "Number"
          data_all = rbind(data_all,data_temp)
        }
        
        
        data_temp = data.frame()
        for (ii in Class3_list_temp){
          subsubset1 = subset(subsubd1, Class3 != ii)
          subsubset2 = subset(subsubd1, Class3 == ii )
          aa = setdiff(subsubset2$reads,subsubset1$reads)
          data_temp = rbind(data_temp,data.frame(HLAname = ii, Num = length(aa),stringsAsFactors = FALSE))
          
        }
        
        if (nrow(data_temp) > 0){
          data_temp$Class = "Class3"
          data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
          data_temp$Genes = geneName
          data_temp$Type = "Unique"
          data_all = rbind(data_all,data_temp)
        }
        
        
        data_temp = data.frame()
        for (ii in Class4_list_temp){
          subsubset2 = subset(subsubd1, Class4 == ii )
          aa = unique(subsubset2$reads )
          data_temp = rbind(data_temp,data.frame(HLAname = ii, Num = length(aa),stringsAsFactors = FALSE))
          
        }
        
        if (nrow(data_temp) > 0){
          data_temp$Class = "Class4"
          data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
          data_temp$Genes = geneName
          data_temp$Type = "Number"
          data_all = rbind(data_all,data_temp)
        }
        
        
        
        data_temp = data.frame()
        for (ii in Class4_list_temp){
          aa = paste(sapply(strsplit(ii, split=':', fixed=TRUE), function(x) (x[1:3])),collapse  = ":")
          subsubset1 = subset(subsubd1, Class4 != ii & Class3 == aa )
          subsubset2 = subset(subsubd1, Class4 == ii )
          aa = setdiff(subsubset2$reads,subsubset1$reads)
          data_temp = rbind(data_temp,data.frame(HLAname = ii, Num = length(aa),stringsAsFactors = FALSE))
          
        }
        
        if (nrow(data_temp) > 0){
          data_temp$Class = "Class4"
          data_temp = data_temp[order(data_temp$Num,decreasing=TRUE),]
          data_temp$Genes = geneName
          data_temp$Type = "Unique"
          data_all = rbind(data_all,data_temp)
        }
        
      }
      
    }
    
  }
}

data_all = subset(data_all,Num>0)

write.csv(data_all,paste0(gsub(".bam.txt|.Nuc.txt","",file_ID),"_HLAlist.csv"),row.names = FALSE)






