#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(mclust)
library(ggplot2)


option_list = list(
   make_option(c("-f", "--file"), type="character", default="pca_output.csv", 
               help="dataset file name [default=%default]", metavar="character"),
   make_option(c("-o", "--out"), type="character", default="ligands_cluster", 
               help="output file name base [default=%default]", metavar="character"),
   make_option(c("-n", "--ncluster"), type="integer", default=11, 
               help="Number of clusters to use for mclust [default=%default]", metavar="integer"),
   make_option(c("-i", "--iter"), type="integer", default=5, 
               help="Number of iterations of mclust for clustering [default=%default]", metavar="integer"),
   make_option(c("-F", "--eff"), type="double", default=0.55, 
               help="Minimum efficicacy of ligand for inclusion in cluster files [default=%default]", metavar="double"),
   make_option(c("-e", "--energy"), type="double", default=-8, 
               help="Maximum binding energy of ligand for inclusion in cluster files [default=%default]", metavar="double"),
   make_option(c("-g", "--graph"), type="character", default="cluster_graph.png", 
               help="Filename of graph of clusters [default=%default]", metavar="character")
); 

opt_parser = OptionParser(usage = "usage: %prog [options]",
                          option_list = option_list,
                          add_help_option = TRUE,
                          prog = NULL,
                          description = "Performs clustering og ligand binding data.",
                          epilogue = "");
opt = parse_args(opt_parser);

zinc<- read.csv(opt$file)

data <- zinc[,c(3, 5, 8, 10, 11, 14:ncol(zinc))]
zinc.pca <- prcomp(data, center = TRUE,scale. = TRUE)
var_pca <- zinc.pca$sdev^2
total_var <- sum(var_pca)
cum_var <- cumsum(var_pca)/total_var
last_ind <- length(cum_var[cum_var<0.9])+1
ind.coord <- as.data.frame(zinc.pca$x)[,c(1:last_ind)]

max_ll = -Inf

subset <- sample(x = 1:nrow(ind.coord), size = 2500,replace = FALSE)
i = 1

while(i <= opt$iter){
   # k <- Mclust(ind.coord, G = opt$ncluster, modelNames = c("VEV","VVV"), initialization = list(hcPairs, subset) )
   # k <- Mclust(ind.coord, G = (opt$ncluster-5):opt$ncluster,  modelNames = c("VEV","VVV"), initialization = list(hcPairs, subset), prior = priorControl(shrinkage = 0.1) )
   k <- Mclust(ind.coord, G = opt$ncluster,  modelName = "VVV", initialization = list(hcRandomPairs(ind.coord)) )
   if(is.null(k)){
      print(paste("Model was NULL. Check subset selection. Retrying with new subset.", sep=""))
   }
   else if(k$loglik > max_ll){
      max_ll = k$loglik 
      clust <- k
      print(paste("Iteration ", i, " of ", opt$iter, ". New best log-likelihood!", sep=""))
      print(summary(clust))
      i=i+1
   }
   else{
      print(paste("Iteration ", i, " of ", opt$iter, ". Log-likelihood of last model (", k$modelName, ") was ", round(k$loglik,3) , ". Best log-likelihood is ", round(max_ll,3) , ".", sep = "" ))
      i = i+1
   }
}


zinc$group <- factor(clust$classification)


ind.coord2 <- ind.coord
ind.coord2$cluster <- zinc$group
ind.coord2$eff <- zinc$Eff
ind.coord2$id <- zinc$ZincID
ind.coord2$energy <- zinc$Energy
ind.coord2$scale <- 0.00001*abs(ind.coord2$energy^4)


ind.coord.filtered <- ind.coord2 %>% filter(eff > opt$eff, energy < opt$energy) %>% droplevels()
zinc.filtered <- zinc %>% filter(Efficiency > opt$eff, Energy < opt$energy) %>% droplevels()

for(i in 1:length(levels(zinc.filtered$group))){
   file.name=paste(opt$out,levels(zinc.filtered$group)[i],".csv", sep = "")
   curr_data = zinc.filtered %>% filter(group == levels(factor(k$classification))[i])
   write.csv(curr_data,file=file.name, quote=FALSE, row.names=FALSE)
}


p<- ind.coord.filtered %>%
   ggplot(aes(x = PC1, y = PC2, color = cluster, label = id)) +
   geom_point(aes(size = scale)) +
   scale_size(name = waiver(), breaks = waiver(), labels = waiver(), limits = NULL, range = c(0.1, 3.5), trans = "identity", guide = "legend") + 
   guides(scale = FALSE, size = FALSE)

ggsave(filename = opt$graph, plot = p,  device = "png", units = "in", width = 10, height = 7)

p.full<- ind.coord2 %>%
   ggplot(aes(x = PC1, y = PC2, color = cluster, label = id)) +
   geom_point(aes(size = scale)) +
   scale_size(name = waiver(), breaks = waiver(), labels = waiver(), limits = NULL, range = c(0.1, 3.5), trans = "identity", guide = "legend") + 
   guides(scale = FALSE, size = FALSE)

ggsave(filename = paste(unlist(strsplit(opt$graph, "\\."))[1], "_full.png", sep=""), plot = p.full,  device = "png", units = "in", width = 10, height = 7)

#library(plotly)
  # ggplotly(p.full, tooltip = c("label","cluster"))
