rm(list = ls())
#setwd("D:\\Projects\\HLA\\test")

url <- "https://raw.githubusercontent.com/ANHIG/IPDKIR/refs/heads/Latest/kir_nuc.fasta"
file_path <- "kir_nuc.fasta"

# Download the file
download.file(url, file_path, mode = "wb")
print(getwd())
if (file.exists(file_path)){
  print("The KIR sequence data has been successfully downloaded.")
  suppressPackageStartupMessages({
	library(Biostrings)
	library(data.table)
	library(dplyr)
	})
  sequences <- readDNAStringSet(file_path)
  
  seq = paste0(">",names(sequences),"\n",sequences[])
  seq = unique(seq)
  
  aa = sapply(strsplit( seq , split=' ', fixed=TRUE), function(x) (x[ ]))
  
  
  data = data.frame(V1 = aa[1,], V2 = aa[2,], V3 = aa[3,], V4 = aa[4,])
  
  data$V1 = gsub(">","",data$V1)
  
  original_strings <-  as.character(data$V2)
  substring_after_asterisk <- substr(original_strings, 
                                     start = regexpr("\\*", original_strings) + 1, 
                                     stop = nchar(original_strings))
  # Count the number of characters in the extracted substring
  data$Number <- nchar(substring_after_asterisk)
  
  
  data_all = data.frame()
  d1 = data%>%
    filter(Number <= 4)
  data_all = rbind(data_all,d1)
  data_all$V5 = data_all$V2
  
  
  d1 = data%>%
    filter(Number == 5)
  modified_strings <- gsub("(\\d{3})(\\d{2})", "\\1:\\2", d1$V2)
  d1$V5 = modified_strings
  data_all = rbind(data_all,d1)
  
  d1 = data%>%
    filter(Number == 6)
  modified_strings <- gsub("(\\d{3})(\\d{2})", "\\1:\\2", d1$V2)
  d1$V5 = modified_strings
  data_all = rbind(data_all,d1)
  
  d1 = data%>%
    filter(Number == 7)
  modified_strings <-  gsub("(\\d{3})(\\d{2})(\\d{2})", "\\1:\\2:\\3", d1$V2)
  d1$V5 = modified_strings
  data_all = rbind(data_all,d1)
  
  d1 = data%>%
    filter(Number == 8)
  modified_strings <-  gsub("(\\d{3})(\\d{2})(\\d{2})", "\\1:\\2:\\3", d1$V2)
  d1$V5 = modified_strings
  data_all = rbind(data_all,d1)
  
  data_all$V5 = paste0("KIR",data_all$V5)
  
  
  
  summary_df <- data_all %>%
    group_by(V4) %>%
    summarise(v1_combined = paste(V5, collapse = ", "), .groups = "drop")%>%
    distinct()
  
  aa = sapply(strsplit( summary_df$v1_combined , split=',', fixed=TRUE), function(x) (x[1 ]))
  data_all = data_all %>%
    filter(V5 %in% aa)
  
  
  d1 = paste0(">",data_all$V5," ",data_all$V3," ",data_all$V4)
  cat(d1, file = paste0("data/",file_path), append = FALSE,sep = "\n") 
  
  d1 = data_all[,c(1,6,2,5)]
  d1$V2 = paste0("KIR",d1$V2)
  fwrite(d1,"data/KIRName.txt",append = FALSE, quote = "auto", sep="\t",row.names = FALSE, col.names = FALSE)
  
  
  d1 = data_all[,c(1,6,3)]
  d1$V3 = as.numeric(as.character(d1$V3))
  fwrite(d1,"data/KIR-length.txt",append = FALSE, quote = "auto", sep="\t",row.names = FALSE, col.names = FALSE)
  
  colnames(summary_df) = c("Sequence","Alleles")
  summary_df$Sequence = gsub("bp\n","",summary_df$Sequence)
  summary_df = summary_df[grepl(",",summary_df$Alleles),]
  write.csv(summary_df, "data/KIR.name.csv",row.names = FALSE)
  file.remove(file_path)
}else{
  print("The KIR sequence data cannot be downloaded at this time.")
  
}





