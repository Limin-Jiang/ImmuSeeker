#   ImmuSeeker
#
#   Authors:    Limin Jiang (lxj423@med.miami.edu)
#   Date:       June 1, 2024
#   Version:    0.1.0

rm(list = ls())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggraph)
  library(igraph)
  library(docopt)
})
#filenames = "P001D00.bam"
#threhold_num = 0
doc <- "
Usage: myscript.R [-i <Imputname>] [-n <Numreads>]
Options:
  -i, --Imputname   The file of list of cancers
  -n, --Numreads   The file of list of cancers
"
options <- docopt(doc)

## parameters
filenames = options$Imputname 
threhold_num = as.numeric(options$Numreads ) 



d1 = fread(paste0(filenames,"_genetype.csv"),header = TRUE,sep = ",")
listss = unique(as.vector(sapply(strsplit(d1$Genetype, split=',', fixed=TRUE), function(x) (x[]))))

d1 = fread(paste0(filenames,"_list.csv"),header = TRUE,sep = ",") %>%
  filter(Type == "Unique")

d1$ID = sapply(strsplit(d1$HLAname, split=':', fixed=TRUE), function(x) (x[1]))
#print(sapply(d1, class))
d1 = d1%>%
  filter((ID %in% listss) & (Num >= threhold_num) )

#print(nrow(d1))
#write.csv(d1,paste0(filenames,"temp.csv"))
getdir = function(subdata_class3,subdata_class4){
  data = data.frame()
  for (ii in 1:nrow(subdata_class3) ){
    HLA_ID = subdata_class3$HLAname[ii]
    subdata = subdata_class4[grepl( paste0("^",  gsub("\\*","\\\\*",HLA_ID),":" ),subdata_class4$HLAname),]
    if (nrow(subdata)){
      data = rbind(data,data.frame(from = paste0(HLA_ID,"(",subdata_class3$Num[ii],")" ) ,to = paste0(subdata$HLAname,"(",subdata$Num,")" )))
    }
    
  }
  return(data)
  
}


getdata = function(d1,gene_ID){
  
  
  subdata = subset(d1, Class == "Class1")
  if (nrow(subdata) > 0){
    subdata_class1 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class1 = NULL
  }
  
  
  subdata = subset(d1, Class == "Class2")
  if (nrow(subdata) > 0){
    subdata_class2 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class2 = NULL
  }
  
  
  subdata = subset(d1, Class == "Class3")
  if (nrow(subdata) > 0){
    subdata_class3 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class3 = NULL
  }
  
  
  
  subdata = subset(d1, Class == "Class4")
  if (nrow(subdata) > 0){
    subdata_class4 = subdata[order(subdata$Num,decreasing=TRUE),]
  }else{
    subdata_class4 = NULL
  }
  
  
  data_all = data.frame()
  
  if (!is.null(subdata_class1)){
    
    data = data.frame(from = gene_ID,to = paste0(subdata_class1$HLAname,"(",subdata_class1$Num,")" ))
    data_all = rbind(data_all,data)
  }
  
  
  
  if (!is.null(subdata_class1) & !is.null(subdata_class2)){
    data_all = rbind(data_all,getdir(subdata_class1,subdata_class2))
    
  }
 
  if (!is.null(subdata_class2) & !is.null(subdata_class3)){
    data_all = rbind(data_all,getdir(subdata_class2,subdata_class3))
  }
  if (!is.null(subdata_class3) & !is.null(subdata_class4)){
    data_all = rbind(data_all,getdir(subdata_class3,subdata_class4))
  }
  return(data_all)
  
}

if (nrow(d1) > 1){
  genes = unique(d1$Genes)
  pdf(paste0(filenames,"_evolution_graphs.pdf"), width = 6, height = 4.5) 
  for (gene_ID in genes){
	  print(gene_ID)
    subd1 = d1%>%
      filter(Genes == gene_ID)
    
    data_all = getdata(subd1,gene_ID)
    
    mygraph <- graph_from_data_frame( data_all )
    
    plot <- ggraph(mygraph, layout = 'auto', circular = FALSE) + 
      geom_edge_diagonal(strength = 1) +
      geom_node_point(color = 'grey', size = 1)+
      geom_node_text(aes(label=name) ,colour = "blue", angle=90 ,size = 1.5, hjust=0.2,vjust=-0.7) +
      theme_void()+
      ggtitle(paste0("HLA ",gene_ID))
    
    print(plot)
  }
}
dev.off()

