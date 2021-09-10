#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(mclust)
library(ggplot2)
library(RColorBrewer)


option_list = list(
   make_option(c("-f", "--file"), type="character", default="bayesMM_ligands.csv", 
               help="dataset file name [default=%default]", metavar="character"),
   make_option(c("-o", "--out"), type="character", default="ligands_cluster", 
               help="output file name base [default=%default]", metavar="character"),
   make_option(c("-F", "--eff"), type="double", default=0.4, 
               help="Minimum efficicacy of ligand for inclusion in cluster files [default=%default]", metavar="double"),
   make_option(c("-e", "--energy"), type="double", default=-8.5, 
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

gg_color_hue <- function(n) {
   hues = seq(15, 375, length = n + 1)
   hues = hcl(h = hues, l = 65, c = 100)[1:n]
   # new_hues = hues[c(1,4,7,2,5,8,3,6)]
   # new_hues
   hues
}

# setwd("/Users/James/Downloads/Screening/test")
zinc<- read.csv(opt$file)

zinc$scale <- 0.00001*abs(zinc$Energy^4)

zinc$group.f <- as.factor(zinc$group)
n.lines <- nlevels(zinc$group.f)
n.colors <- ceiling(n.lines/3)
cols = gg_color_hue(n.colors)
colors <- rep(cols, times = n.lines%/%n.colors + 1)
shapes <- rep(c(19, 17, 15), each = n.colors)

zinc$PC1.s <- sign(zinc$PC1)*sqrt(1-exp(-zinc$PC1^2/9))
zinc$PC2.s <- sign(zinc$PC2)*sqrt(1-exp(-zinc$PC2^2/9))
zinc.filtered <- zinc %>% group_by(group.f) %>% slice_min(order_by = Energy, n =200)


minx = min(zinc$PC1.s)
maxx = max(zinc$PC1.s)

miny = min(zinc$PC2.s)
maxy = max(zinc$PC2.s)


for(i in 1:max(zinc$group)){
   file.name=paste(opt$out,i,".csv", sep = "")
   curr_data = zinc %>% filter(group.f == i)
   write.csv(curr_data,file=file.name, quote=FALSE, row.names=FALSE)

   p.clust <- zinc %>% filter(group == i) %>%
      ggplot(aes(x = PC1.s, y = PC2.s, color = group.f, label = ZincID)) +
      geom_point(aes(size = scale, shape = group.f)) +
      scale_size(name = waiver(), breaks = waiver(), labels = waiver(), limits = NULL, range = c(0.1, 2.5), trans = "identity", guide = "legend") + 
      guides(scale = FALSE, size = FALSE) + 
      scale_color_manual(values=colors[i]) +
      scale_shape_manual(values=shapes[i]) +
      xlim(minx, maxx) + 
      ylim(miny, maxy)

   ggsave(filename = paste("cluster_", i, ".png", sep=""), plot = p.clust,  device = "png", units = "in", width = 10, height = 7)


}


p <- zinc.filtered %>%
   ggplot(aes(x = PC1.s, y = PC2.s, color = group.f, label = ZincID)) +
   geom_point(aes(size = scale, shape = group.f)) +
   scale_size(name = waiver(), breaks = waiver(), labels = waiver(), limits = NULL, range = c(0.1, 2.5), trans = "identity", guide = "legend") + 
   guides(scale = FALSE, size = FALSE) + 
   scale_color_manual(values=colors) +
   scale_shape_manual(values=shapes) 
# p

# ggplotly(p, tooltip = c("label","group"))

ggsave(filename = opt$graph, plot = p,  device = "png", units = "in", width = 10, height = 7)

p.full<- zinc %>% 
   ggplot(aes(x = PC1.s, y = PC2.s, color = group.f, label = ZincID)) +
   geom_point(aes(size = scale, shape = group.f)) +
   scale_size(name = waiver(), breaks = waiver(), labels = waiver(), limits = NULL, range = c(0.1, 2.5), trans = "identity", guide = "legend") + 
   guides(scale = FALSE, size = FALSE) + 
   scale_color_manual(values=colors) +
   scale_shape_manual(values=shapes) 

ggsave(filename = paste(unlist(strsplit(opt$graph, "\\."))[1], "_full.png", sep=""), plot = p.full,  device = "png", units = "in", width = 10, height = 7)



# library(plotly)
# ggplotly(p.full, tooltip = c("label","group"))
