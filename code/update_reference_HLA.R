rm(list = ls())
library(data.table)
#setwd("D:\\Projects\\HLA\\test")

url <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_nuc.fasta"
destfile <- "hla_nuc.fasta"

# Download the file
download.file(url, destfile, mode = "wb")

if (file.exists(destfile)){
  print("The HLA sequence data has been successfully downloaded.")
  
  library(Biostrings)
  sequences <- readDNAStringSet(destfile)
  
  seq = paste0(">",names(sequences),"\n",sequences[])
  seq = unique(seq)
  
  aa = sapply(strsplit( seq , split=' ', fixed=TRUE), function(x) (x[ ]))
  d1 = paste0(">HLA-",aa[2,]," ",aa[3,]," ",aa[4,])
  
  cat(d1, file = paste0("./data/",destfile), sep = "\n")
  
  
  d1 = data.frame(V1 = gsub(">","",aa[1,]), V2 = aa[2,], V3 = as.numeric(aa[3,])  )
  fwrite(d1, file = "./data/HLA-length.txt", append = FALSE, quote = "auto",
         sep="\t",row.names = FALSE, col.names = FALSE)
  
  
  aa = sapply(strsplit( unique(names(sequences)) , split=' ', fixed=TRUE), function(x) (x[ ]))
  d1 = data.frame(V1 =aa[1,],V2 =aa[2,] )
  fwrite(d1, file = "./data/HLA_name.txt", append = FALSE, quote = "auto",
         sep="\t",row.names = FALSE, col.names = FALSE)
  
  file.remove(destfile)
  
}else{
  print("The HLA sequence data cannot be downloaded at this time.")
  
}













