###################### This script shows how to make plots from tidygraph objects downloaded from SpeSpeNet.https://doi.org/10.1016/j.watres.2022.118767

library(ggplot2)
library(pals)
library(tidygraph)
library(ggraph)
library(ggdark)
library(igraph)

### Define plotting themes:

themes <- list("Dark" = dark_theme_grey(),
               "Classic" = theme_bw(),
               "Minimal" = theme_minimal())

lineCol <- list("Dark" = "white",
                "Classic" = "black",
                "Minimal" = "black")

use_theme <- theme(axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   legend.title = element_text(size = 28),
                   legend.text = element_text(size = 30),
                   legend.position = c(0, 0),
                   legend.justification = c("left", "bottom"),
                   legend.background = element_rect(fill='transparent'),
                   legend.margin = margin(6, 6, 6, 6))


########################## If the network is overlayed with a continuous environmental variable and only positive correlations are used:


### Load in downloaded tidygraph object:

tidy.net <- readRDS(file = "Tidygraph_object_continuous_env.rds")

#### Optional: Delete isolated nodes

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

#### Rotate the network until it doesnt overlap with the legend
set.seed(3)

plot.env <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = envcor),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = bquote(O[2]*" correlation")) +
  guides(fill = guide_colorbar(barheight = 10, barwidth = 5)) +
  themes[[3]]+
  use_theme+
  theme(legend.text = element_text(size =20),legend.position = c(0,0.02))

png(filename = "net_env_continuous.png",width = 500,height = 500)
plot.env
dev.off()

########################## If the network is overlayed with a continuous environmental variable and postive and negative correlations are used:

tidy.net <- readRDS(file = "Tidygraph_object_continuous_env_bothcor.rds")

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

set.seed(5)

plot.env <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    aes(edge_color = sign)) +
  geom_point(aes(x,y,size = nodeSize,fill = envcor),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  scale_edge_color_manual(values = c("#FFA500","#0047AB")) + 
  themes[[3]] +
  labs(fill = bquote(O[2]*" correlation")) +
  guides(fill = guide_colorbar(barheight = 10, barwidth = 5),edge_color = "none") +
  use_theme+
  theme(legend.title.align = 0.5,legend.text = element_text(size = 18),legend.position = c(0,0.02))

png(filename = "net_env_continuous_allcor.png",width = 500,height = 500)
plot.env
dev.off()


########################## If the network is overlayed with a discrete environmental variable and postive correlations are used:

tidy.net <- readRDS(file = "Tidygraph_object_discrete.rds")

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

### Manually define color palette equal to number of categories in the variable
col.palette <- c("#A6CEE3", "#1F78B4", "#33A02C", "#FFFF00", "#FDBF6F", "#FB9A99", "#E31A1C")

set.seed(3)

plot.max <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = envcor),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_manual(values = col.palette) +
  themes[["Classic"]] +
  guides(fill = guide_legend(title = "Maximum mean\nrelative abundance in:",override.aes = list(size = 7),title.position = "top")) +
  use_theme+
  theme(legend.title.align = 0.5,legend.text = element_text(size = 18))

png(filename = "net_env_discrete.png",width = 500,height = 500)
plot.max
dev.off()

########################## If the network is overlayed with a discrete environmental variable and postive correlations are used:

tidy.net <- readRDS(file = "Tidygraph_object_discrete_env_bothcor.rds")

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

### Manually define color palette equal to number of categories in the variable
col.palette <- c("#A6CEE3", "#1F78B4", "#33A02C", "#FFFF00", "#FDBF6F", "#FB9A99", "#E31A1C")

set.seed(5)

plot.max <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    aes(edge_color = sign)) +
  geom_point(aes(x,y,size = nodeSize,fill = envcor),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_manual(values = col.palette) +
  scale_edge_color_manual(values = c("#FFA500","#0047AB")) + 
  themes[["Classic"]] +
  guides(fill = guide_legend(title = "Maximum mean\nrelative abundance in:",override.aes = list(size = 7),title.position = "top"),edge_color = "none") +
  use_theme+
  theme(legend.title.align = 0.5,legend.text = element_text(size = 18))

png(filename = "net_env_discrete_allcor.png",width = 500,height = 500)
plot.max
dev.off()

########################## If the network is overlayed with taxonomy and positive correlations are used:

tidy.net <- readRDS(file = "Tidygraph_object_taxonomy.rds")

node.data <- tidy.net%>%activate(nodes)%>%data.frame()
tax.data <- node.data$taxo

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

set.seed(3)

plot.tax <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = taxo),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend(title = "Class", override.aes = list(size = 7),ncol = 1)) +
  scale_fill_manual(values = as.vector(c(polychrome(length(unique(tax.data)))[-c(1,2)],"#E4E1E3","#5A5156")), na.value = "grey") +
  themes[["Classic"]] +
  use_theme +
  theme(legend.text = element_text(size = 15),legend.title = element_text(size = 20))

png(filename = "net_taxonomy.png",width = 500,height = 500)
plot.tax
dev.off()

########################## If the network is overlayed with taxonomy and positive and negative correlations are used:

tidy.net <- readRDS(file = "Tidygraph_object_taxonomy_bothcor.rds")

node.data <- tidy.net%>%activate(nodes)%>%data.frame()
tax.data <- node.data$taxo

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

set.seed(3)

plot.tax <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    aes(edge_color = sign)) +
  geom_point(aes(x,y,size = nodeSize,fill = taxo),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend(title = "Class", override.aes = list(size = 7),ncol = 1),edge_color = "none") +
  scale_fill_manual(values = as.vector(c(polychrome(length(unique(tax.data)))[-c(1,2)],"#E4E1E3","#5A5156")), na.value = "grey") +
  scale_edge_color_manual(values = c("#FFA500","#0047AB")) + 
  themes[["Classic"]] +
  use_theme +
  theme(legend.text = element_text(size = 15),legend.title = element_text(size = 20))

png(filename = "net_taxonomy_allcor.png",width = 500,height = 500)
plot.tax
dev.off()

########################## If the network is overlayed with kmeans cluster and positive correlations are used:

tidy.net <- readRDS(file = "Tidygraph_object_cluster.rds")

node.data <- tidy.net%>%activate(nodes)%>%data.frame()
clus.data <- node.data$cluster

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

set.seed(3)

plot.clust <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = as.factor(cluster)),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend("Cluster", override.aes = list(size = 15))) +
  scale_fill_manual(values = as.vector(polychrome(length(unique(clus.data))+4)[-c(1,2,4,5)])) +
  themes[["Classic"]] +
  use_theme

png(filename = "net_cluster.png",width = 500,height = 500)
plot.clust
dev.off()

########################## If the network is overlayed with kmeans cluster and positive and negative correlations are used:

tidy.net <- readRDS(file = "Tidygraph_object_cluster_bothcor.rds")

node.data <- tidy.net%>%activate(nodes)%>%data.frame()
clus.data <- node.data$cluster

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

set.seed(3)

plot.clust <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    aes(edge_color = sign)) +
  geom_point(aes(x,y,size = nodeSize,fill = as.factor(cluster)),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_edge_color_manual(values = c("#FFA500","#0047AB")) + 
  guides(fill = guide_legend("Cluster", override.aes = list(size = 15)),edge_color = "none") +
  scale_fill_manual(values = as.vector(polychrome(length(unique(clus.data))+4)[-c(1,2,4,5)])) +
  themes[["Classic"]] +
  use_theme

png(filename = "net_cluster_allcor.png",width = 500,height = 500)
plot.clust
dev.off()






