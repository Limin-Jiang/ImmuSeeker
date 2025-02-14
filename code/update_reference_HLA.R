rm(list = ls())
#setwd("D:\\Projects\\HLA\\test")

url <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_nuc.fasta"
destfile <- "hla_nuc.fasta"

# Download the file
download.file(url, destfile, mode = "wb")

if (file.exists(destfile)){
  print("The HLA sequence data has been successfully downloaded.")
  
  suppressPackageStartupMessages({
	library(Biostrings)
	library(data.table)
	library(dplyr)
	})
  sequences <- readDNAStringSet(destfile)
  
  seq = paste0(">",names(sequences),"\n",sequences[])
  seq = unique(seq)
  
  aa = sapply(strsplit( seq , split=' ', fixed=TRUE), function(x) (x[ ]))
  
  
  data_all = data.frame(V1 = gsub(">","", aa[1,]), V2 = aa[2,], V3 = aa[3,], V4 = aa[4,] )
  
  summary_df <- data_all %>%
    group_by(V4) %>%
    summarise(v1_combined = paste(V2, collapse = ", "), .groups = "drop")%>%
    distinct()
  
  aa = sapply(strsplit( summary_df$v1_combined , split=',', fixed=TRUE), function(x) (x[1 ]))
  data_all = data_all %>%
    filter(V2 %in% aa)
  
  
  colnames(summary_df) = c("Sequence","Alleles")
  summary_df$Sequence = gsub("bp\n","",summary_df$Sequence)
  summary_df = summary_df[grepl(",",summary_df$Alleles),]
  write.csv(summary_df, "data/HLA.name.csv",row.names = FALSE)

  
  d1 = paste0(">HLA-",data_all$V2," ",data_all$V3," ",data_all$V4)
  cat(d1, file = paste0("data/",destfile), sep = "\n")
  
  
  
  d1 = data.frame(V1 = gsub(">","",data_all$V1), V2 = data_all$V2, V3 = as.numeric(data_all$V3)  )
  fwrite(d1, file = "data/HLA-length.txt", append = FALSE, quote = "auto",
         sep="\t",row.names = FALSE, col.names = FALSE)
  
  
  aa = sapply(strsplit( unique(names(sequences)) , split=' ', fixed=TRUE), function(x) (x[ ]))
  d1 = data.frame(V1 =aa[1,],V2 =aa[2,] )
  fwrite(d1, file = "data/HLA_name.txt", append = FALSE, quote = "auto",
         sep="\t",row.names = FALSE, col.names = FALSE)
  
  file.remove(destfile)
  
}else{
  print("The HLA sequence data cannot be downloaded at this time.")
  
}













