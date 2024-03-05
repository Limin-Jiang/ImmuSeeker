#   HLApipeline
#
#   Authors:    Limin Jiang (lxj423@med.miami.edu)
#   Date:       2024-3-1
#   Version:    0.2.0
# Default values


rm(list = ls())
library(tidyverse)
library(phyloseq) #otu_table
library(rstatix) #t_test
library(ggprism) #add_pvalue
library(ggplot2)

###################Alpha diversity Example##################################

setwd("C:\\Users\\lxj423\\OneDrive - University of Miami\\Desktop\\pipline")
loaded_data  = readRDS("diversity_analysis_data.rds")
Label = loaded_data$Label
data_end = loaded_data$data_end

OTU = otu_table(data_end, taxa_are_rows = TRUE)
richness <- estimate_richness(OTU,measures = c("Shannon", "Simpson", "InvSimpson"))


for (ii in 1:ncol(richness)){
  data_end = data.frame(Label = Label, Value = richness[,ii],stringsAsFactors = FALSE)
  alpha_ID = colnames(richness)[ii]
  stat.test <- data_end %>%
    t_test(Value ~ Label) %>%
    add_significance()%>% 
    add_y_position() 
  
  if ( nrow(stat.test) > 0){
    each.vs.ref <- tibble::tribble(
      ~group1, ~group2, ~p.adj, ~y.position,
      stat.test$group1[1], stat.test$group2[1],stat.test$p[1],stat.test$y.position[1]*(1.01),
      stat.test$group1[2], stat.test$group2[2],stat.test$p[2],stat.test$y.position[2]*(1.03),
      stat.test$group1[3], stat.test$group2[3],stat.test$p[3],stat.test$y.position[3]*(1.05)
    )
    each.vs.ref = each.vs.ref %>%
      filter(p.adj < 0.05)
    each.vs.ref$p.adj = signif(each.vs.ref$p.adj,1)
    each.vs.ref$p.adj[each.vs.ref$p.adj >= 0.0001] = paste("p =",each.vs.ref$p.adj[each.vs.ref$p.adj >= 0.0001],sep = " ")
    each.vs.ref$p.adj[as.numeric(each.vs.ref$p.adj)  < 0.0001] = "p < 0.0001"
    each.vs.ref$p.adj[each.vs.ref$p.adj == "p = 1e-04"] = "p = 0.0001"
    
    colorlist = data.frame(Control = "#6E8B3D", MS = "#CD3700", NMO = "#009ACD")
    ggplot(data_end, aes(x = Label, y = Value)) +  
      geom_boxplot(aes(fill = Label, color= Label),alpha=0.5)+ 
      geom_point(aes(fill = Label), position = position_jitter(width = 0.2),shape = 25,size = 2)+
      add_pvalue(each.vs.ref) +
      theme_classic()+
      theme(
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", size = 0.01, color = "Gainsboro"),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 12,color = "black",face = "bold"),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title = element_text(size = 12,face = "bold"),
      ) +
      scale_color_manual(values = colorlist)+
      scale_fill_manual(values =colorlist)+
      xlab("") +
      ylab(paste0(alpha_ID," Index"))
    path <- paste0(alpha_ID,"_","alpha", ".jpeg")
    ggsave(path, width =3, height = 3)
    
  }
  
  
}

###################beta diversity Example##################################
rm(list = ls())
library(vegan)
library(ggplot2)
library(tidyverse)
setwd("C:\\Users\\lxj423\\OneDrive - University of Miami\\Desktop\\pipline")
loaded_data  = readRDS("diversity_analysis_data.rds")
Label = loaded_data$Label
data_end = loaded_data$data_end %>%
  t()%>%
  as.data.frame()
  


methods = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger")
for (lag_ID in methods[1]){ ######This is a example#########
  bray_curtis_matrix <- vegdist(data_end, method = lag_ID)
  nmds_result <- metaMDS(bray_curtis_matrix)
  
  Value_points =  nmds_result$points %>%
    as.data.frame()
  Value_points$Label = Label
  
  colorlist = data.frame(Control = "#6E8B3D", MS = "#CD3700", NMO = "#009ACD")
  ggplot(data=Value_points, aes(x = MDS1 , y = MDS2,color = Label )) + 
    geom_point(size=3,alpha = 0.7) +
    stat_ellipse(geom = "polygon",
                 aes(fill = Label), 
                 alpha = 0.25,
                 linetype = 2) +
    theme_bw()+
    scale_color_manual(values = colorlist)+
    scale_fill_manual(values =colorlist)+
    theme(axis.title = element_text(size = 14,color = "black",face = "bold"),
          axis.text.x = element_text(size = 12,color = "black", angle = 0, hjust=0.5,vjust=0.5),
          axis.text.y = element_text(size = 12,color = "black", hjust = 0.2),
          panel.grid = element_blank(),
          legend.position="top",
          legend.title = element_text( size = 14),
          legend.text = element_text( size = 14,color = "black",face = "bold"),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))+
    guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))+
    labs(y = "MDS2", x = "MDS1")
  
  
  
  path <- paste0(lag_ID,"_BETA_diversity.jpeg")
  ggsave(path, width =4, height = 4)
}













