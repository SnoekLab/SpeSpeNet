.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6")

library(tidyverse)
library(tools)
library(igraph)
library(pals)
library(reshape2)
library(tidygraph)
library(scales)
library(ggraph)
library(ggdark)

relab <- read.delim("NW_data/Brenzinger/otu.new.txt", header = T, sep = "\t",row.names = 1)
tax <- read.delim("NW_data/Brenzinger/tax.new.txt", header = T, sep = "\t",row.names = 1)
tax[is.na(tax)] <- "Unknown"
tax[tax == "NA"] <- "Unknown"
colnames(tax) <- toTitleCase(colnames(tax))

relab <- aggregate(relab,tax[,tolower(colnames(tax))!="species"],sum)
tax <- relab %>% dplyr::select(where(is.character))
relab <- relab %>%  dplyr::select(where(is.numeric))

relab.rel <- sweep(relab, 2, colSums(relab), FUN="/") * 100

mmo <- NULL
mmo$mean <- apply(relab.rel, 1, mean, na.rm = T)
mmo$max <- apply(relab.rel, 1, max, na.rm = T)
mmo$occ <- apply(relab.rel > 0, 1, sum, na.rm = T)

relab <- relab[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 10),]
relab.rel <- relab.rel[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 10),]
tax <- tax[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 10),]

set.seed(333)
relab[relab==0] <- runif(sum(relab==0),0.1,1)

row.name <- rownames(relab)

relab <- apply(relab,2,function(x)log(x/mean(x)))

rownames(relab) <- row.name

use_relab <- relab

env <- read.delim("NW_data/Brenzinger/meta.new.txt", header = T, sep = "\t",row.names = 1)
env$Soil_type[which(env$Soil_type=="sandy")] <- "Sand"
env$Soil_type[which(env$Soil_type=="clay")] <- "Clay"

env_num <- env %>%
  select(where(is.numeric))
env_num <- env_num[,which(apply(env_num,2,function(x)length(unique(x))>1))]
env_fac <- env %>% 
  select(where(is.factor)|where(is.character))
env_fac <- env_fac[which(apply(env_fac,2,function(x) length(unique(x))>1))]
env_sub <- env_fac[which(apply(env_fac,2,function(x) max(table(x),na.rm = T)>7))]

cor_mat <- cor(base::t(use_relab), method = "spearman")

relab.long <- melt(data.frame(use_relab,otu = rownames(use_relab)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["amendment"]][match(relab.long$variable,rownames(env))])
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)
colnames(relab.mean) <- c("otu","sub_cat","value")
max.cat <- sapply(unique(relab.mean$otu),function(x) relab.mean[relab.mean$otu==x,"sub_cat"][which.max(relab.mean[relab.mean$otu==x,"value"])])
env.cor.amendment <- max.cat[match(rownames(use_relab),names(max.cat))]
env.cor.amendment[env.cor.amendment=="control65"]  <- "Control"

relab.long <- melt(data.frame(use_relab,otu = rownames(use_relab)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["Soil_type"]][match(relab.long$variable,rownames(env))])
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)
colnames(relab.mean) <- c("otu","sub_cat","value")
max.cat <- sapply(unique(relab.mean$otu),function(x) relab.mean[relab.mean$otu==x,"sub_cat"][which.max(relab.mean[relab.mean$otu==x,"value"])])
env.cor.soiltype <- max.cat[match(rownames(use_relab),names(max.cat))]
env.cor.ph <- as.numeric(cor(env_num[["pH"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.om <- as.numeric(cor(env_num[["OM"]],base::t(use_relab), use = "pairwise",method = "pearson"))

diag(cor_mat) <- NA

km.mat <- cor_mat
diag(km.mat) <- 1

set.seed(1)
km.clus <- kmeans(km.mat, 2, nstart = 25)

cor.thr <- 0.79

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
  slice_max(sum, n = 9) %>%
  pull("Order")

# Change taxonomy to "Other" for entries outside of the top 11
final.tax <- case_when(tax.mrg[["Phylum"]] %in% tax.top ~ tax.mrg[["Phylum"]], TRUE ~ "Other")
final.tax[final.tax=="Unclass_Bacteria"] <- "Unknown"

# Change taxonomy to "Other" for entries outside of the top 11
final.tax.order <- case_when(tax.mrg.order[["Order"]] %in% tax.top.order ~ tax.mrg.order[["Order"]], TRUE ~ "Other")


tidy.net <- otu_net %>%
  as_tbl_graph()
# Node size is proportional to the square root of the abundance of the corresponding OTU
meanSize <- sqrt(rowMeans(relab.rel)) %>%
  rescale(to = c(1,12))

tidy.net %<>%
  activate(nodes) %>%
  mutate(nodeSize = meanSize,
         cluster = km.clus$cluster,
         env_cor_amendment = env.cor.amendment,
         env_cor_soiltype = env.cor.soiltype,
         env_cor_ph = env.cor.ph,
         env_cor_om = env.cor.om,
         taxo = factor(final.tax, levels = c(sort(unique(final.tax[!(final.tax%in%c("Other","Unknown"))])),"Other","Unknown")),
         counts = rowMeans(use_relab),
         order = final.tax.order)

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

set.seed(14)

plot.max.amendment <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_amendment),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_brewer(palette = "Paired", na.value = "grey") +
  themes[["Classic"]] +
  guides(fill = guide_legend(title = "Preferred\namendment",override.aes = list(size = 7),title.position = "top")) +
  use_theme+
  theme(legend.title.align = 0,legend.text = element_text(size = 22),legend.title = element_text(size = 21),legend.position = c(0,-0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_full_amendment.png",width = 500,height = 500)
plot.max.amendment
dev.off()

set.seed(14)

plot.max.st <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_soiltype),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_manual(values = c("red","blue"), na.value = "grey") +
  themes[["Classic"]] +
  guides(fill = guide_legend(title = "Preferred\nsoil type",override.aes = list(size = 7),title.position = "top")) +
  use_theme+
  theme(legend.title.align = 0,legend.text = element_text(size = 22),legend.title = element_text(size = 21),legend.position = c(0,0))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_full_soiltype.png",width = 500,height = 500)
plot.max.st
dev.off()

######## Figure S4 ##########

set.seed(14)

s5.ph <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_ph),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\n with pH") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 5)) +
  use_theme+
  theme(legend.text = element_text(size =18),legend.title = element_text(size = 21),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_full_ph.png",width = 500,height = 500)
s5.ph
dev.off()

set.seed(14)

s5.om <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_om),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\n with OM") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 5)) +
  use_theme+
  theme(legend.text = element_text(size =18),legend.title = element_text(size = 21),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_full_om.png",width = 500,height = 500)
s5.om
dev.off()

s5 <- plot_grid(s5.ph,s5.om,plot.max.st,labels = c("A","B","C"),label_size = 35,label_colour = "black",nrow = 2)

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/s5.png",width = 1000,height = 1000)
s5
dev.off()

############ Only clay samples ##########

relab <- read.delim("NW_data/Brenzinger/otu.new.txt", header = T, sep = "\t",row.names = 1)
tax <- read.delim("NW_data/Brenzinger/tax.new.txt", header = T, sep = "\t",row.names = 1)
tax[is.na(tax)] <- "Unknown"
tax[tax == "NA"] <- "Unknown"
colnames(tax) <- toTitleCase(colnames(tax))

relab <- aggregate(relab,tax[,tolower(colnames(tax))!="species"],sum)
tax <- relab %>% select(where(is.character))
relab <- relab %>%  select(where(is.numeric))

env <- read.delim("NW_data/Brenzinger/meta.new.txt", header = T, sep = "\t",row.names = 1)
env$Soil_type[which(env$Soil_type=="sandy")] <- "Sand"
env$Soil_type[which(env$Soil_type=="clay")] <- "Clay"
relab <- relab[,env$Soil_type=="Clay"]
env <- env[env$Soil_type=="Clay",]

env_num <- env %>%
  select(where(is.numeric))
env_num <- env_num[,which(apply(env_num,2,function(x)length(unique(x))>1))]
env_fac <- env %>% 
  select(where(is.factor)|where(is.character))
env_fac <- env_fac[which(apply(env_fac,2,function(x) length(unique(x))>1))]
env_sub <- env_fac[which(apply(env_fac,2,function(x) max(table(x),na.rm = T)>7))]

relab.rel <- sweep(relab, 2, colSums(relab), FUN="/") * 100

mmo <- NULL
mmo$mean <- apply(relab.rel, 1, mean, na.rm = T)
mmo$max <- apply(relab.rel, 1, max, na.rm = T)
mmo$occ <- apply(relab.rel > 0, 1, sum, na.rm = T)

relab <- relab[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 10),]
relab.rel <- relab.rel[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 10),]
tax <- tax[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 10),]

set.seed(333)
relab[relab==0] <- runif(sum(relab==0),0.1,1)

row.name <- rownames(relab)

relab <- apply(relab,2,function(x)log(x/mean(x)))

rownames(relab) <- row.name

use_relab <- relab

cor_mat <- cor(base::t(use_relab), method = "spearman")

env.cor.nirk <- as.numeric(cor(env_num[["nirK"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.plant_heights <- as.numeric(cor(env_num[["plant_heights"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.nifh <- as.numeric(cor(env_num[["nifH"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.noszII <- as.numeric(cor(env_num[["nosZII"]],base::t(use_relab), use = "pairwise",method = "pearson"))

relab.long <- melt(data.frame(use_relab,otu = rownames(use_relab)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["amendment"]][match(relab.long$variable,rownames(env))])
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)
colnames(relab.mean) <- c("otu","sub_cat","value")
max.cat <- sapply(unique(relab.mean$otu),function(x) relab.mean[relab.mean$otu==x,"sub_cat"][which.max(relab.mean[relab.mean$otu==x,"value"])])
env.cor.amendment <- max.cat[match(rownames(use_relab),names(max.cat))]
env.cor.amendment[env.cor.amendment=="control65"]  <- "Control"

diag(cor_mat) <- NA

km.mat <- cor_mat
diag(km.mat) <- 1

set.seed(1)
km.clus <- kmeans(km.mat, 2, nstart = 25)

cor.thr <- 0.66

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
means <- data.frame(name = rownames(relab.rel),
                    meanAbundance = rowMeans(relab.rel))
tax.mrg <- merge(means, tax.flt, by = "name", sort = FALSE)

# Merge graph nodes, mean abundances and taxonomy
nodes.order <- data.frame(name = V(otu_net)$name)
tax.flt.order <- tax.sel.order %>% filter(name %in% nodes$name)
means.order <- data.frame(name = rownames(relab.rel),
                          meanAbundance = rowMeans(relab.rel))
tax.mrg.order <- merge(means, tax.flt.order, by = "name", sort = FALSE)

# Create list of top 11 most abundant taxonomies
tax.top <- tax.mrg %>%
  group_by(.data[["Phylum"]]) %>%
  summarise(sum = sum(meanAbundance)) %>%
  slice_max(sum, n = 8) %>%
  pull("Phylum")

tax.top.order <- tax.mrg.order %>%
  group_by(.data[["Order"]]) %>%
  summarise(sum = sum(meanAbundance)) %>%
  slice_max(sum, n = 9) %>%
  pull("Order")

# Change taxonomy to "Other" for entries outside of the top 11
final.tax <- case_when(tax.mrg[["Phylum"]] %in% tax.top ~ tax.mrg[["Phylum"]], TRUE ~ "Other")
final.tax[final.tax=="Unclass_Bacteria"] <- "Unknown"

# Change taxonomy to "Other" for entries outside of the top 11
final.tax.order <- case_when(tax.mrg.order[["Order"]] %in% tax.top.order ~ tax.mrg.order[["Order"]], TRUE ~ "Other")

tidy.net <- otu_net %>%
  as_tbl_graph()
# Node size is proportional to the square root of the abundance of the corresponding OTU
meanSize <- sqrt(rowMeans(relab.rel)) %>%
  rescale(to = c(1,12))

tidy.net %<>%
  activate(nodes) %>%
  mutate(nodeSize = meanSize,
         cluster = km.clus$cluster,
         env_cor_amendment = env.cor.amendment,
         env_cor_nirk = env.cor.nirk,
         env_cor_nifh = env.cor.nifh,
         env_cor_plant_heights = env.cor.plant_heights,
         env_cor_nosII = env.cor.noszII,
         taxo = factor(final.tax, levels = c(sort(unique(final.tax[!(final.tax%in%c("Other","Unknown"))])),"Other","Unknown")),
         counts = rowMeans(relab.rel),
         order = final.tax.order)

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

set.seed(9)

plot.max.clay.amendment <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_amendment),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_brewer(palette = "Paired", na.value = "grey") +
  themes[["Classic"]] +
  guides(fill = guide_legend(title = "Preferred\namendment",override.aes = list(size = 7),title.position = "top")) +
  use_theme+
  theme(legend.title.align = 0,legend.text = element_text(size = 21),legend.title = element_text(size = 21),legend.position = c(0,-0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_clay_amendment.png",width = 500,height = 500)
plot.max.clay.amendment
dev.off()

set.seed(9)

plot.clay.nirk <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_nirk),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\nwith nirK") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 4)) +
  use_theme+
  theme(legend.text = element_text(size =18),legend.title = element_text(size = 21),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_clay_nirk.png",width = 500,height = 500)
plot.clay.nirk
dev.off()

set.seed(9)

plot.tax <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = taxo),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 6),ncol = 1)) +
  scale_fill_manual(values = as.vector(c(polychrome(length(unique(final.tax))+1)[-c(1,2)],"#E4E1E3")), na.value = "grey") +
  themes[["Classic"]] +
  use_theme +
  theme(legend.text = element_text(size = 19),legend.position = c(-0.01,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_phylum.png",width = 500,height = 500)
plot.tax
dev.off()

set.seed(9)

plot.clay.plant.heights <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_plant_heights),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation with\nplant height") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 4)) +
  use_theme+
  theme(legend.text = element_text(size =18),legend.title = element_text(size = 21),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_clay_plant_heights.png",width = 500,height = 500)
plot.clay.plant.heights
dev.off()

set.seed(9)

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
  theme(text = element_text(size = 25),legend.text = element_text(size = 21),legend.title = element_text(size = 21))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_k2.png",width = 500,height = 500)
plot.clust
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
#net_data <- rbind(data.frame(net_data),data.frame("1","Thaumarchaeota",0,0))

barplot.phylum <- ggplot(net_data, aes(x = cluster, y = Percent_cluster, fill = taxo)) +
  geom_bar(stat = "identity", color = "black") +
  themes[["Classic"]] +
  scale_fill_manual(values = as.vector(c(polychrome(length(unique(final.tax))+1)[-c(1,2)],"#E4E1E3"))) +
  labs(x = "Cluster", y = "Percentage",fill = "Phylum") +
  theme(text = element_text(size = 30),legend.position = "none")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/barplot_clay_phylum.png",width = 300,height = 500)
barplot.phylum
dev.off()

cor.frame <- tbl.net %>%
  as_tibble()
cor.frame$cluster <- factor(cor.frame$cluster)
colnames(cor.frame)[colnames(cor.frame)=="taxo"] <- "Taxonomy"
cor.plot<-ggplot(cor.frame,aes(x = cluster, y = env_cor_nirk)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = Taxonomy),size = 1.5) +
  scale_color_manual(values = as.vector(c(polychrome(length(unique(final.tax)))[-c(1,2)],"#E4E1E3","#5A5156"))) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth = .5) +
  themes[["Classic"]] +
  labs(x = "Cluster", y = "nirK correlation") +
  theme(text=element_text(size=30),legend.position = "none")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/corplot_phylum.png",width = 300,height = 500)
cor.plot
dev.off()

library(cowplot)

#################### Make figure 1:

bar.cor <- plot_grid(barplot.phylum,cor.plot,labels = c("F","G"),label_size = 35,label_colour = "black")

fig_1 <-plot_grid(plot.max.amendment,plot.max.st,plot.max.clay.amendment,plot.clust,plot.tax,bar.cor,plot.clay.nirk,plot.clay.plant.heights,
                      align = "none",nrow = 3, labels = c("A","B","C","D","E","","H","I"),label_size = 35,label_colour = "black",label_x = 0.02)

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/Figure_1.png",width = 1600,height = 1500)
fig_1
dev.off()

#################### Make figure S3

set.seed(9)

plot.clay.nifh <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_nifh),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\n with nifH") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 4)) +
  use_theme+
  theme(legend.text = element_text(size =18),legend.title = element_text(size = 21),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_nifh.png",width = 500,height = 500)
plot.clay.nifh
dev.off()

set.seed(9)

plot.clay.nosZII <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_nosII),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\n with nosZII") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 4)) +
  use_theme+
  theme(legend.text = element_text(size =18),legend.title = element_text(size = 21),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_nosZII.png",width = 500,height = 500)
plot.clay.nosZII
dev.off()

fig.s3 <- plot_grid(plot.clay.nifh,plot.clay.nosZII,labels = c("A","B"),label_size = 35,label_colour = "black")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/Figure_s3.png",width = 1000,height = 500)
fig.s3
dev.off()

#################### Make Figure S1

set.seed(9)

plot.clay.nosZII.default <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Dark"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_nosII),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Dark"]] +
  labs(fill = "Correlation\n with nosZII") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 4,ticks.colour = NA)) +
  use_theme+
  theme(legend.text = element_text(size =16),legend.title = element_text(size = 18),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_nosZII_default.png",width = 500,height = 500)
plot.clay.nosZII.default
dev.off()

set.seed(16)

plot.clay.nosZII.altered.iso <- ggraph(tbl.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.7,
    alpha = 0.6,
    strength = 0.8,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_nosII),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\n with nosZII") +
  guides(fill = guide_colorbar(barheight = 8, barwidth = 4,ticks.colour = NA)) +
  use_theme+
  theme(legend.text = element_text(size =16),legend.title = element_text(size = 18),legend.position = "none")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/net_nosZII_altered_iso.png",width = 500,height = 500)
plot.clay.nosZII.altered.iso
dev.off()

fig.s1 <- plot_grid(plot_grid(plot.clay.nosZII.default,labels = c("A"),label_size = 35,label_colour = "white"),
                        plot_grid(plot.clay.nosZII.altered.iso,labels = c("B"),label_size = 35,label_colour = "black"))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Brenzinger_CLS/Figure_S1.png",width = 1000,height = 500)
fig.s1
dev.off()






