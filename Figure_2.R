.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6")

library(tidyverse)
library(tools)
library(igraph)
library(pals)
library(reshape2)

relab <- read.delim("NW_data/Meisner/Meisner.relab.txt", header = T, sep = "\t",row.names = 1)
tax <- read.delim("NW_data/Meisner/altered_Meisner.tax.txt", header = T, sep = "\t",row.names = 1)

tax[is.na(tax)] <- "Unknown"
tax[tax == "NA"] <- "Unknown"
colnames(tax) <- toTitleCase(colnames(tax))

env <- read.delim("NW_data/Meisner/altered_Meisner.meta.txt", header = T, sep = "\t",row.names = 1)

env_num <- env %>%
  select(where(is.numeric))
env_num <- env_num[,which(apply(env_num,2,function(x)length(unique(x))>1))]
env_fac <- env %>% 
  select(where(is.factor)|where(is.character))
env_fac <- env_fac[which(apply(env_fac,2,function(x) length(unique(x))>1))]
env_sub <- env_fac[which(apply(env_fac,2,function(x) max(table(x),na.rm = T)>7))]
mmo <- NULL
mmo$mean <- apply(relab, 1, mean, na.rm = T)
mmo$max <- apply(relab, 1, max, na.rm = T)
mmo$occ <- apply(relab > 0, 1, sum, na.rm = T)

use_relab <- relab[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 5),]
tax <- tax[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 5),]

cor_mat <- cor(base::t(use_relab), method = "spearman")

env.cor <- as.numeric(cor(env_num[["Legacy.num"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.norg <- as.numeric(cor(env_num[["Ninorg"]],base::t(use_relab), use = "pairwise",method = "pearson"))

diag(cor_mat) <- NA

km.mat <- cor_mat
diag(km.mat) <- 1

set.seed(1)

km.clus <- kmeans(km.mat, 5, nstart = 25)

cor.thr <- 0.5

otu_net <- graph.adjacency(cor_mat, weighted = TRUE, mode = "upper")
otu_net <- delete.edges(otu_net, which(E(otu_net)$weight < cor.thr | E(otu_net)$weight %in% NA))
otu_net <- igraph::simplify(otu_net, remove.multiple = T, remove.loops = T)

tax.sel <- tax %>% select("Phylum")
tax.sel <- tibble::rownames_to_column(tax.sel, "name")

tax.sel.order <- tax %>% select("Order")
tax.sel.order <- tibble::rownames_to_column(tax.sel.order, "name")

# Merge graph nodes, mean abundances and taxonomy
nodes <- data.frame(name = V(otu_net)$name)
tax.flt <- tax.sel %>% filter(name %in% nodes$name)
means <- data.frame(name = rownames(use_relab),
                    meanAbundance = rowMeans(use_relab))
tax.mrg <- merge(means, tax.flt, by = "name", sort = FALSE)

# Merge graph nodes, mean abundances and taxonomy
nodes.order <- data.frame(name = V(otu_net)$name)
tax.flt.order <- tax.sel.order %>% filter(name %in% nodes$name)
means.order <- data.frame(name = rownames(use_relab),
                    meanAbundance = rowMeans(use_relab))
tax.mrg.order <- merge(means, tax.flt.order, by = "name", sort = FALSE)

# Create list of top 11 most abundant taxonomies
tax.top <- tax.mrg %>%
  group_by(.data[["Phylum"]]) %>%
  summarise(sum = sum(meanAbundance)) %>%
  slice_max(sum, n = 9) %>%
  pull("Phylum")

tax.top.order <- tax.mrg.order %>%
  group_by(.data[["Order"]]) %>%
  summarise(sum = sum(meanAbundance)) %>%
  slice_max(sum, n = 13) %>%
  pull("Order")

tax.top.order <- c(tax.top.order,"Sphingomonadales")

# Change taxonomy to "Other" for entries outside of the top 11
final.tax <- case_when(tax.mrg[["Phylum"]] %in% tax.top ~ tax.mrg[["Phylum"]], TRUE ~ "Other")
final.tax[final.tax=="Unclass_Bacteria"] <- "Unknown"

# Change taxonomy to "Other" for entries outside of the top 11
final.tax.order <- case_when(tax.mrg.order[["Order"]] %in% tax.top.order ~ tax.mrg.order[["Order"]], TRUE ~ "Other")


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
           envcornorg = env.cor.norg,
           taxo = factor(final.tax, levels = c(sort(unique(final.tax[!(final.tax%in%c("Other","Unknown"))])),"Other","Unknown")),
           counts = rowMeans(use_relab),
           order = final.tax.order)
} else{
  tidy.net %<>%
    activate(nodes) %>%
    mutate(nodeSize = meanSize,
           cluster = km.clus$cluster,
           taxo = final.tax,
           counts = rowMeans(use_relab))
}

tbl.net <- tidy.net

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)
if(length(rem_nodes)!=0){
  final.tax <- final.tax[-c(rem_nodes)]
}

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

set.seed(29)

### Km plot
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
  use_theme+
  theme(text = element_text(size = 25))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Meisner/net_k5.png",width = 500,height = 500)
plot.clust
dev.off()

set.seed(29)

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
  labs(fill = "Moisture\ncorrelation") +
  guides(fill = guide_colorbar(barheight = 10, barwidth = 5)) +
  use_theme+
  theme(legend.text = element_text(size =20))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Meisner/net_legacy_num.png",width = 500,height = 500)
plot.env
dev.off()

set.seed(29)

plot.norg<- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = envcornorg),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Inorganic N\ncorrelation") +
  guides(fill = guide_colorbar(barheight = 10, barwidth = 5)) +
  use_theme+
  theme(legend.text = element_text(size =20))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Meisner/net_norg.png",width = 500,height = 500)
plot.norg
dev.off()

set.seed(29)

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
  theme(legend.text = element_text(size = 23))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Meisner/net_phylum.png",width = 500,height = 500)
plot.tax
dev.off()

net_data <- tbl.net %>%
  activate(nodes) %>%
  as_tibble() %>%
  group_by(cluster, order) %>%
  tally(counts) %>%
  mutate(freq = n / sum(n) * 100) %>%
  ungroup()
net_data$cluster <- factor(net_data$cluster)
colnames(net_data)[which(colnames(net_data)=="n")]<-"Percent_total"
colnames(net_data)[which(colnames(net_data)=="freq")]<-"Percent_cluster"

net_data$order[net_data$order=="Unclass_Bacteria"] <- "Unknown"
net_data$order[net_data$order=="Unclass_Acidobacteria_Gp6"] <- "Acidobacteria_Gp6 (unclassified)"
net_data$order[net_data$order=="Unclass_Actinobacteria"] <- "Actinobacteria (unclassified)"
net_data$order[net_data$order=="Unclass_Bacteroidetes"] <- "Bacteroidetes (unclassified)"
net_data$order[net_data$order=="Unclass_Betaproteobacteria"] <- "Betaproteobacteria (unclassified)"
net_data$order[net_data$order=="Unclass_Gammaproteobacteria"] <- "Gammaproteobacteria (unclassified)"
net_data$order <- factor(net_data$order, levels = c(sort(unique(net_data$order[!(net_data$order%in%c("Other","Unknown"))])),"Other","Unknown"))

barplot.order <- ggplot(net_data, aes(x = cluster, y = Percent_total, fill = order)) +
  geom_bar(stat = "identity", color = "black") +
  themes[["Classic"]] +
  scale_fill_manual(values = as.vector(c(polychrome(length(unique(final.tax.order)))[-c(1,2)],"#E4E1E3","#5A5156"))) +
  labs(x = "Cluster", y = "Percentage",fill = "Phylum") +
  theme(text = element_text(size = 30),legend.position = "none")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Meisner/barplot_order.png",width = 500,height = 500)
barplot.order
dev.off()

cor.frame <- tbl.net %>%
  as_tibble()
cor.frame$cluster <- factor(cor.frame$cluster)
colnames(cor.frame)[colnames(cor.frame)=="order"] <- "Taxonomy"

cor.frame$Taxonomy[cor.frame$Taxonomy=="Unclass_Bacteria"] <- "Unknown"
cor.frame$Taxonomy[cor.frame$Taxonomy=="Unclass_Acidobacteria_Gp6"] <- "Acidobacteria_Gp6 (unclassified)"
cor.frame$Taxonomy[cor.frame$Taxonomy=="Unclass_Actinobacteria"] <- "Actinobacteria (unclassified)"
cor.frame$Taxonomy[cor.frame$Taxonomy=="Unclass_Bacteroidetes"] <- "Bacteroidetes (unclassified)"
cor.frame$Taxonomy[cor.frame$Taxonomy=="Unclass_Betaproteobacteria"] <- "Betaproteobacteria (unclassified)"
cor.frame$Taxonomy[cor.frame$Taxonomy=="Unclass_Gammaproteobacteria"] <- "Gammaproteobacteria (unclassified)"
cor.frame$Taxonomy <- factor(cor.frame$Taxonomy, levels = c(sort(unique(cor.frame$Taxonomy[!(cor.frame$Taxonomy%in%c("Other","Unknown"))])),"Other","Unknown"))

cor.plot<-ggplot(cor.frame,aes(x = cluster, y = envcornorg)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = Taxonomy),size = 1.5) +
  scale_color_manual(values = as.vector(c(polychrome(length(unique(final.tax.order)))[-c(1,2)],"#E4E1E3","#5A5156"))) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth = .5) +
  themes[["Classic"]] +
  labs(x = "Cluster", y = "Inorganic N correlation") +
  theme(text=element_text(size=30),legend.position = "none")

cor.plot

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Meisner/corplot_order.png",width = 300,height = 500)
cor.plot
dev.off()

bar.cor.legend <- get_legend(ggplot(net_data, aes(x = cluster, y = Percent_total, fill = order)) +
                               geom_bar(stat = "identity", color = "black") +
                               themes[["Classic"]] +
                               scale_fill_manual(values = as.vector(c(polychrome(length(unique(final.tax.order)))[-c(1,2)],"#E4E1E3","#5A5156"))) +
                               labs(x = "Cluster", y = "Percentage",fill = "") +
                               theme(text = element_text(size = 35)))

library(cowplot)
library(lemon)

bar.cor <- plot_grid(barplot.order,cor.plot,labels = c("E","F"),label_size = 35,label_colour = "black")

plot.comb <-plot_grid(plot.env,plot.tax,plot.norg,plot.clust,bar.cor,bar.cor.legend,align = "none",nrow = 2,
                      labels = c("A","B","C","D","",""),label_size = 35,label_colour = "black",label_x = 0)

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Meisner/comb_plot.png",width = 1600,height = 1000)
plot.comb
dev.off()



