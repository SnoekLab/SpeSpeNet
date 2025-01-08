###################### This scripts generates figure 2 of the paper titled "SpeSpeNet: An interactive and user-friendly tool 
###################### to create and explore microbial correlation networks." The data used is from Hauptfeld et al 2022 
###################### https://doi.org/10.1016/j.watres.2022.118767

library(tidyverse)
library(tools)
library(igraph)
library(pals)
library(reshape2)
library(tidygraph)
library(scales)
library(ggraph)
library(ggdark)
library(cowplot)
library(SpiecEasi)

############### Load and process data

relab <- read.delim("NW_data/Griftpark/grift.otu.counts.txt", header = T, sep = "\t",row.names = 1)
tax <- read.delim("NW_data/Griftpark/grift.tax.counts.txt", header = T, sep = "\t",row.names = 1)

############### Aggregate at genus level

relab <- aggregate(relab,tax[,tolower(colnames(tax))!="species"],sum)
tax <- relab %>% select(where(is.character))
relab <- relab %>%  select(where(is.numeric))

############### Set unknown taxa as unknown

tax[is.na(tax)] <- "Unknown"
tax[tax == "NA"] <- "Unknown"
colnames(tax) <- toTitleCase(colnames(tax))

############### Load in environmental data

env <- read.delim("NW_data/Griftpark/fixed.grift.meta.txt", header = T, sep = "\t",row.names = 1)

env$sampling_point <- c("Park side: Park area","Park side: Park area","Park side: Collection basin","Park side: Collection basin","Park side: Collection basin",
                        "Park side: Collection basin",rep("Park side: Pipeline",5),rep("Plant side: Influent buffer",4),rep("Plant side: Anti-bulking reactor",4),rep("Plant side: Over-flow",4),
                        rep("Plant side: Sludge-thickener",2))

env$sampling_point <- factor(env$sampling_point,levels = c("Park side: Park area","Park side: Collection basin","Park side: Pipeline","Plant side: Influent buffer",
                                                           "Plant side: Anti-bulking reactor","Plant side: Over-flow","Plant side: Sludge-thickener"))

############### Calculate relative abundances

relab.rel <- sweep(relab, 2, colSums(relab), FUN="/") * 100

############### Filter data

mmo <- NULL
mmo$mean <- apply(relab.rel, 1, mean, na.rm = T)
mmo$max <- apply(relab.rel, 1, max, na.rm = T)
mmo$occ <- apply(relab.rel > 0, 1, sum, na.rm = T)

use_relab <- relab[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 5),]
tax <- tax[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 5),]
relab.rel <- relab.rel[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 5),]

############### Add pseudo-counts

set.seed(333)
use_relab[use_relab==0] <- runif(sum(use_relab==0),0.1,1)

############### Calculate correlation matrix

set.seed(333)
cor_mat <- t(sparcc(t(use_relab)))[[2]]
rownames(cor_mat) <- rownames(use_relab)
colnames(cor_mat) <- rownames(use_relab)

############### Calculate correlation with numerical variable

cls.norm <- use_relab
row.name <- rownames(use_relab)
cls.norm <- apply(cls.norm,2,function(x)log(x/mean(x)))
rownames(cls.norm) <- row.name

env.cor <- as.numeric(cor(env[["O2.levels"]],base::t(cls.norm), use = "pairwise",method = "pearson"))

############### Calculate in which environment taxa have highest mean abundance

relab.long <- melt(data.frame(use_relab,otu = rownames(use_relab)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["sampling_point"]][match(relab.long$variable,rownames(env))])
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)
colnames(relab.mean) <- c("otu","sub_cat","value")
max.cat <- sapply(unique(relab.mean$otu),function(x) relab.mean[relab.mean$otu==x,"sub_cat"][which.max(relab.mean[relab.mean$otu==x,"value"])])
env.cor.max <- max.cat[match(rownames(use_relab),names(max.cat))]

############### Do k-means clustering on correlation matrix

diag(cor_mat) <- NA

km.mat <- cor_mat
diag(km.mat) <- 1

set.seed(1)
km.clus <- kmeans(km.mat, 2, nstart = 25)

############### Make network

cor.thr <- 0.5

otu_net <- graph.adjacency(cor_mat, weighted = TRUE, mode = "upper")
otu_net <- delete.edges(otu_net, which(E(otu_net)$weight < cor.thr | E(otu_net)$weight %in% NA))
otu_net <- igraph::simplify(otu_net, remove.multiple = T, remove.loops = T)

############### Get order rank taxonomy for each node

tax.sel <- tax %>% select("Order")
tax.sel <- tibble::rownames_to_column(tax.sel, "name")

nodes <- data.frame(name = V(otu_net)$name)
tax.flt <- tax.sel %>% filter(name %in% nodes$name)
means <- data.frame(name = rownames(relab.rel),
                    meanAbundance = rowMeans(relab.rel))
tax.mrg <- merge(means, tax.flt, by = "name", sort = FALSE)

###############  Create list of top 9 most abundant taxonomies

tax.top <- tax.mrg %>%
  group_by(.data[["Order"]]) %>%
  summarise(sum = sum(meanAbundance)) %>%
  slice_max(sum, n = 9) %>%
  pull("Order")

############### Change taxonomy to "Other" for entries outside of the top 9
final.tax <- case_when(tax.mrg[["Order"]] %in% tax.top ~ tax.mrg[["Order"]], TRUE ~ "Other")
final.tax[final.tax=="unknown"] <- "Unknown"

############### Make tidynet object

tidy.net <- otu_net %>%
  as_tbl_graph()
# Node size is proportional to the square root of the abundance of the corresponding OTU
meanSize <- sqrt(rowMeans(use_relab)) %>%
  rescale(to = c(1,12))
if(!is.null(env.cor)){
  tidy.net %<>%
    activate(nodes) %>%
    mutate(nodeSize = meanSize,
           cluster = km.clus$cluster,
           envcor = env.cor,
           taxo = factor(final.tax, levels = c(sort(unique(final.tax[!(final.tax%in%c("Other","Unknown"))])),"Other","Unknown")),
           counts = rowMeans(use_relab),
           envmax = env.cor.max)
} else{
  tidy.net %<>%
    activate(nodes) %>%
    mutate(nodeSize = meanSize,
           cluster = km.clus$cluster,
           taxo = final.tax,
           counts = rowMeans(use_relab))
}

tbl.net <- tidy.net

############### Delete isolated nodes

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)
if(length(rem_nodes)!=0){
  final.tax <- final.tax[-c(rem_nodes)]
}

############### Make network plots

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
  scale_fill_manual(values = as.vector(polychrome(length(unique(km.clus$cluster))+3)[-c(1,2,3)])) +
  themes[["Classic"]] +
  use_theme

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Grift_sparcc/net_k2.png",width = 500,height = 500)
plot.clust
dev.off()

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
  use_theme+
  theme(legend.text = element_text(size =20),legend.position = c(0,0.02))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Grift_sparcc/net_O2.png",width = 500,height = 500)
plot.env
dev.off()

set.seed(3)

plot.tax <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = taxo),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 7),ncol = 1)) +
  scale_fill_manual(values = as.vector(c(polychrome(length(unique(final.tax)))[-c(1,2)],"#E4E1E3","#5A5156")), na.value = "grey") +
  themes[["Classic"]] +
  use_theme +
  theme(legend.text = element_text(size = 18))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Grift_sparcc/net_order.png",width = 500,height = 500)
plot.tax
dev.off()

set.seed(3)

col.palette <- c("#A6CEE3", "#1F78B4", "#33A02C", "#FFFF00", "#FDBF6F", "#FB9A99", "#E31A1C")

plot.max <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = envmax),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_manual(values = col.palette) +
  themes[["Classic"]] +
  guides(fill = guide_legend(title = "Maximum mean\nrelative abundance in:",override.aes = list(size = 7),title.position = "top")) +
  use_theme+
  theme(legend.title.align = 0.5,legend.text = element_text(size = 18))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Grift_sparcc/net_stage.png",width = 500,height = 500)
plot.max
dev.off()

net_data <- tbl.net %>%
  activate(nodes) %>%
  as_tibble() %>%
  group_by(cluster, taxo) %>%
  tally(counts) %>%
  mutate(freq = n / sum(n) * 100) %>%
  ungroup()
net_data$cluster <- factor(net_data$cluster)
colnames(net_data)[which(colnames(net_data)=="n")]<-"Percent_total"
colnames(net_data)[which(colnames(net_data)=="freq")]<-"Percent_cluster"

barplot.order <- ggplot(net_data, aes(x = cluster, y = Percent_cluster, fill = taxo)) +
  geom_bar(stat = "identity", color = "black") +
  themes[["Classic"]] +
  scale_fill_manual(values = as.vector(c(polychrome(length(unique(final.tax)))[-c(1,2)],"#E4E1E3","#5A5156"))) +
  labs(x = "Cluster", y = "Percentage",fill = "Order") +
  theme(text = element_text(size = 30),legend.position = "none")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Grift_sparcc/barplot_order.png",width = 300,height = 500)
barplot.order
dev.off()

############### Make barplot and correlation plot per cluster

cor.frame <- tbl.net %>%
  as_tibble()
cor.frame$cluster <- factor(cor.frame$cluster)
colnames(cor.frame)[colnames(cor.frame)=="taxo"] <- "Taxonomy"
cor.plot<-ggplot(cor.frame,aes(x = cluster, y = envcor)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = Taxonomy),size = 1.5) +
  scale_color_manual(values = as.vector(c(polychrome(length(unique(final.tax)))[-c(1,2)],"#E4E1E3","#5A5156"))) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth = .5) +
  themes[["Classic"]] +
  labs(x = "Cluster", y = bquote(O[2]*" correlation")) +
  theme(text=element_text(size=30),legend.position = "none")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Grift_sparcc/corplot_order.png",width = 300,height = 500)
cor.plot
dev.off()

############### Paste all subplots together to generate figure 1.

bar.cor <- plot_grid(barplot.order,cor.plot,labels = c("E","F"),label_size = 35,label_colour = "black")

fig_1 <-plot_grid(plot.tax,plot.env,plot.max,plot.clust,bar.cor,align = "none",nrow = 2,
                      labels = c("A","B","C","D",""),label_size = 35,label_colour = "black",label_x = 0.5)

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Grift_sparcc/Figure_1.png",width = 1600,height = 1000)
fig_1
dev.off()

