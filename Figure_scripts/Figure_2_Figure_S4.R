###################### This scripts generates Figure 2 and Figure S2 of the paper titled "SpeSpeNet: An interactive and  
###################### user-friendly tool to create and explore microbial correlation networks." The data used is from 
###################### Deutschmann et al 2023 https://doi.org/10.1186/s40168-023-01523-z

library(tidyverse)
library(tools)
library(igraph)
library(pals)
library(reshape2)
library(tidygraph)
library(scales)
library(ggraph)
library(ggdark)
library(RColorBrewer)
library(cowplot)

############### Load and process data

relab <- read.delim("NW_data/TemporalNetworkBBMO/asv.txt", header = T, sep = "\t",row.names = 1)
tax <- read.delim("NW_data/TemporalNetworkBBMO/tax.txt", header = T, sep = "\t",row.names = 1)
relab <- sweep(relab, 2, colSums(relab), FUN="/") * 100

tax[is.na(tax)] <- "Unknown"
tax[tax == "NA"] <- "Unknown"
colnames(tax) <- toTitleCase(colnames(tax))

relab <- aggregate(relab,tax[,tolower(colnames(tax))!="species"],sum)
tax <- relab %>% select(where(is.character))
relab <- relab %>%  select(where(is.numeric))

env <- read.delim("NW_data/TemporalNetworkBBMO/env.txt", header = T, sep = "\t",row.names = 1)

env$Seasons <- factor(env$Seasons,levels = c("Winter","Spring","Summer","Autumn"))

############### Filter data

mmo <- NULL
mmo$mean <- apply(relab, 1, mean, na.rm = T)
mmo$max <- apply(relab, 1, max, na.rm = T)
mmo$occ <- apply(relab > 0, 1, sum, na.rm = T)

use_relab <- relab[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 5),]
tax <- tax[which(mmo$mean >= 0 & mmo$max >= 0.1 & mmo$occ >= 5),]

############### Calculate correlation matrix

cor_mat <- cor(base::t(use_relab), method = "spearman")

############### Calculate correlation with numerical variable

env.cor.temp <- as.numeric(cor(env[["Temp"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.chlo <- as.numeric(cor(env[["CHL_total"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.day.length <- as.numeric(cor(env[["Day_length_Hours_light"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.secchi <- as.numeric(cor(env[["SECCHI"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.NO3 <- as.numeric(cor(env[["NO3"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.PO4 <- as.numeric(cor(env[["PO4"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.NO2 <- as.numeric(cor(env[["NO2"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.NH4 <- as.numeric(cor(env[["NH4"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.SI <- as.numeric(cor(env[["SI"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.year <- as.numeric(cor(env[["year"]],base::t(use_relab), use = "pairwise",method = "pearson"))
env.cor.SAL_CTD <- as.numeric(cor(env[["SAL_CTD"]],base::t(use_relab), use = "pairwise",method = "pearson"))

############### Calculate in which environments taxa have highest mean abundance

relab.long <- melt(data.frame(use_relab,otu = rownames(use_relab)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["Seasons"]][match(relab.long$variable,rownames(env))])
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)
colnames(relab.mean) <- c("otu","sub_cat","value")
max.cat <- sapply(unique(relab.mean$otu),function(x) relab.mean[relab.mean$otu==x,"sub_cat"][which.max(relab.mean[relab.mean$otu==x,"value"])])
env.cor.seasons <- max.cat[match(rownames(use_relab),names(max.cat))]

relab.long <- melt(data.frame(use_relab,otu = rownames(use_relab)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["month"]][match(relab.long$variable,rownames(env))])
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)
colnames(relab.mean) <- c("otu","sub_cat","value")
max.cat <- sapply(unique(relab.mean$otu),function(x) relab.mean[relab.mean$otu==x,"sub_cat"][which.max(relab.mean[relab.mean$otu==x,"value"])])
env.cor.month <- max.cat[match(rownames(use_relab),names(max.cat))]

############### Make network

diag(cor_mat) <- NA

cor.thr <- 0.4
otu_net <- graph.adjacency(cor_mat, weighted = TRUE, mode = "upper")
otu_net <- delete.edges(otu_net, which(E(otu_net)$weight < cor.thr | E(otu_net)$weight %in% NA))
otu_net <- igraph::simplify(otu_net, remove.multiple = T, remove.loops = T)

############### Get class rank taxonomy for each node

tax.sel <- tax %>% select("Class")
tax.sel <- tibble::rownames_to_column(tax.sel, "name")

nodes <- data.frame(name = V(otu_net)$name)
tax.flt <- tax.sel %>% filter(name %in% nodes$name)
means <- data.frame(name = rownames(use_relab),
                    meanAbundance = rowMeans(use_relab))
tax.mrg <- merge(means, tax.flt, by = "name", sort = FALSE)

tax.top <- tax.mrg %>%
  group_by(.data[["Class"]]) %>%
  summarise(sum = sum(meanAbundance)) %>%
  slice_max(sum, n = 9) %>%
  pull("Class")

# Change taxonomy to "Other" for entries outside of the top 9
final.tax <- case_when(tax.mrg[["Class"]] %in% tax.top ~ tax.mrg[["Class"]], TRUE ~ "Other")

############### Get kingdom rank taxonomy for each node

tax.kingdom <- tax%>%select("Kingdom")
tax.kingdom <- tibble::rownames_to_column(tax.kingdom,"name")

nodes <- data.frame(name = V(otu_net)$name)
tax.flt <- tax.kingdom %>% filter(name %in% nodes$name)
means <- data.frame(name = rownames(use_relab),
                    meanAbundance = rowMeans(use_relab))
tax.mrg <- merge(means, tax.flt, by = "name", sort = FALSE)

tax.top <- tax.mrg %>%
  group_by(.data[["Kingdom"]]) %>%
  summarise(sum = sum(meanAbundance)) %>%
  slice_max(sum, n = 9) %>%
  pull("Kingdom")

final.tax.king <- case_when(tax.mrg[["Kingdom"]] %in% tax.top ~ tax.mrg[["Kingdom"]], TRUE ~ "Other")

############### Make tidynet object

tidy.net <- otu_net %>%
  as_tbl_graph()
# Node size is proportional to the square root of the abundance of the corresponding OTU
meanSize <- sqrt(rowMeans(use_relab)) %>%
  rescale(to = c(1,12))

tidy.net %<>%
  activate(nodes) %>%
  mutate(nodeSize = meanSize,
         env_cor_month = env.cor.month,
         env_cor_seasons = env.cor.seasons,
         env_cor_chlo = env.cor.chlo,
         env_cor_day_length = env.cor.day.length,
         env_cor_temp = env.cor.temp,
         env_cor_year = env.cor.year,
         env_cor_secchi = env.cor.secchi,
         env_cor_no3 = env.cor.NO3,
         env_cor_po4 = env.cor.PO4,
         env_cor_nh4 = env.cor.NH4,
         env_cor_no2 = env.cor.NO2,
         env_cor_si = env.cor.SI,
         env_cor_sal_CTD = env.cor.SAL_CTD,
         taxo = factor(final.tax, levels = c(sort(unique(final.tax[!(final.tax%in%c("Other","Unknown"))])),"Other","Unknown")),
         taxo_king = factor(final.tax.king, levels = c(sort(unique(final.tax.king[!(final.tax.king%in%c("Other","Unknown"))])),"Other","Unknown")),
         counts = rowMeans(use_relab))

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

set.seed(6)

plot.full.secchi <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_secchi),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\nwith\nSecchi depth") +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_secchi.png",width = 500,height = 500)
plot.full.secchi
dev.off()

set.seed(6)

plot.full.no3 <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_no3),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = bquote(NO[3]^"-"*" correlation")) +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_no3.png",width = 500,height = 500)
plot.full.no3
dev.off()

set.seed(6)

plot.full.chlo <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_chlo),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\nwith\nchlorophyll a") +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_chlorophyll.png",width = 500,height = 500)
plot.full.chlo
dev.off()


set.seed(6)

plot.full.temp <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_temp),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\nwith\ntemperature") +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_temperature.png",width = 500,height = 500)
plot.full.temp
dev.off()

set.seed(6)

plot.full.daylength <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_day_length),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\nwith\nday length") +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_daylength.png",width = 500,height = 500)
plot.full.daylength
dev.off()

set.seed(6)

plot.full.year <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_year),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\nwith\nsampling year") +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_year.png",width = 500,height = 500)
plot.full.year
dev.off()

set.seed(6)

plot.full.PO4 <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_po4),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = bquote(PO[4]^"3-"*" correlation")) +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_PO4.png",width = 500,height = 500)
plot.full.PO4
dev.off()

set.seed(6)

plot.full.NO2 <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_no2),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = bquote(NO[2]*" correlation")) +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_NO2.png",width = 500,height = 500)
plot.full.NO2
dev.off()

set.seed(6)

plot.full.NH4 <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_nh4),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = bquote(NH[4]^"+"*" correlation")) +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_NH4.png",width = 500,height = 500)
plot.full.NH4
dev.off()

set.seed(6)

plot.full.sal <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_sal_CTD),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = "Correlation\nwith\nsalinity") +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_salinity.png",width = 500,height = 500)
plot.full.sal
dev.off()

set.seed(6)

plot.full.si <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_si),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
  themes[["Classic"]] +
  labs(fill = bquote(SiO[2]*" correlation")) +
  guides(fill = guide_colorbar(barheight = 9, barwidth = 4.5)) +
  use_theme+
  theme(legend.text = element_text(size = 18,vjust = 0.7),legend.title = element_text(size = 20),legend.position = c(0,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_si.png",width = 500,height = 500)
plot.full.si
dev.off()

set.seed(6)

plot.tax.suphy <- ggraph(tidy.net, layout = 'fr') +
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
  theme(legend.text = element_text(size = 21))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_superphylum.png",width = 500,height = 500)
plot.tax.suphy
dev.off()

set.seed(6)

plot.tax.king <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = taxo_king),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 7),ncol = 1)) +
  scale_fill_brewer(palette = "Dark2")+
  #scale_fill_manual(values = as.vector(c(polychrome(length(unique(final.tax.king))+2)[-c(1,2)])), na.value = "grey") +
  themes[["Classic"]] +
  use_theme +
  theme(legend.text = element_text(size = 25))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_kingdom.png",width = 500,height = 500)
plot.tax.king
dev.off()

set.seed(6)

plot.max.seasons <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = env_cor_seasons),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  scale_fill_manual(values = c("#A6CEE3","#33A02C","#B2DF8A","#1F78B4")) +
  themes[["Classic"]] +
  guides(fill = guide_legend(title = "Preferred\nseason:",override.aes = list(size = 7),title.position = "top")) +
  use_theme+
  theme(legend.title.align = 0,legend.text = element_text(size = 22),legend.title = element_text(size = 25))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/net_full_seasons.png",width = 500,height = 500)
plot.max.seasons
dev.off()

############### Make separate barplots of distribution per month and per season for the prokaryotic and eukaryotic communities.

tax.16s <- tax[tax$Kingdom!="Eukaryota",]
tax.18s <- tax[tax$Kingdom=="Eukaryota",]

tax.18s[tax.18s=="Arthracanthida-Symphyacanthida"] <- "Acantharea S-A"

use_relab.16s <- use_relab[tax$Kingdom!="Eukaryota",]
use_relab.18s <- use_relab[tax$Kingdom=="Eukaryota",]

############### Rescale counts to 100 in both tables

use_relab.16s <- sweep(use_relab.16s, 2, colSums(use_relab.16s), FUN="/") * 100
use_relab.18s <- sweep(use_relab.18s, 2, colSums(use_relab.18s), FUN="/") * 100

############### Barplots 16S:

taxi <- tax.16s[["Order"]]
names(taxi) <- rownames(tax.16s)
relab.long <- melt(data.frame(use_relab.16s,otu = rownames(use_relab.16s)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["month"]][match(relab.long$variable,rownames(env))])
taxi <- taxi[match(relab.long$otu,names(taxi))]
agg.tax <- aggregate(relab.long$value,by = list(taxi),sum)
most.abu <- if(nrow(agg.tax)>12){
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)][1:12]
} else{
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)]
}
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)

colnames(relab.mean) <- c("otu","sub_cat","value")
use.tax <- tax.16s[match(relab.mean$otu,rownames(tax.16s)),]
relab.mean <- data.frame(relab.mean,use.tax)
relab.mean[["Order"]][!relab.mean[["Order"]]%in%c(most.abu,"unknown")] <- "Other"
relab.mean.agg  <- aggregate(relab.mean[,c("value")],by = list(relab.mean$sub_cat,relab.mean[["Order"]]),sum)
colnames(relab.mean.agg) <- c("Category","Taxonomy","Percent_cluster")
relab.mean.agg$Category <- factor(relab.mean.agg$Category, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
relab.mean.agg$Taxonomy <- gsub("_X","",relab.mean.agg$Taxonomy,fixed = T)
relab.mean.agg$Taxonomy <- gsub("Gammaproteobacteria_","",relab.mean.agg$Taxonomy,fixed = T)
relab.mean.agg$Taxonomy <- factor(relab.mean.agg$Taxonomy, levels = c(sort(unique(relab.mean.agg$Taxonomy[!(relab.mean.agg$Taxonomy%in%c("Other","Unknown"))])),"Other","Unknown"))

month.barplot.16s <- ggplot(relab.mean.agg, aes(x = Category, y = Percent_cluster, fill = Taxonomy)) +
  geom_bar(stat = "identity", color = "black") +
  themes[["Classic"]] +
  scale_fill_manual(values = as.vector(polychrome(15)[-c(1,2)]), na.value = "grey", name = NULL) +
  labs(x = "", y = "Order") +
  ylim(0,100)+
  theme(text = element_text(size = 22),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size =20))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/month_barplot_16s.png",width = 500,height = 500)
month.barplot.16s
dev.off() 

taxi <- tax.16s[["Order"]]
names(taxi) <- rownames(tax.16s)
relab.long <- melt(data.frame(use_relab.16s,otu = rownames(use_relab.16s)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["Seasons"]][match(relab.long$variable,rownames(env))])
taxi <- taxi[match(relab.long$otu,names(taxi))]
agg.tax <- aggregate(relab.long$value,by = list(taxi),sum)
most.abu <- if(nrow(agg.tax)>12){
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)][1:12]
} else{
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)]
}
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)

colnames(relab.mean) <- c("otu","sub_cat","value")
use.tax <- tax.16s[match(relab.mean$otu,rownames(tax.16s)),]
relab.mean <- data.frame(relab.mean,use.tax)
relab.mean[["Order"]][!relab.mean[["Order"]]%in%c(most.abu,"unknown")] <- "Other"
relab.mean.agg  <- aggregate(relab.mean[,c("value")],by = list(relab.mean$sub_cat,relab.mean[["Order"]]),sum)
colnames(relab.mean.agg) <- c("Category","Taxonomy","Percent_cluster")
relab.mean.agg$Category <- factor(relab.mean.agg$Category, levels = c("Winter","Spring","Summer","Autumn"))
relab.mean.agg$Taxonomy <- factor(relab.mean.agg$Taxonomy, levels = c(sort(unique(relab.mean.agg$Taxonomy[!(relab.mean.agg$Taxonomy%in%c("Other","Unknown"))])),"Other","Unknown"))

season.barplot.16s <- ggplot(relab.mean.agg, aes(x = Category, y = Percent_cluster, fill = Taxonomy)) +
  geom_bar(stat = "identity", color = "black") +
  themes[["Classic"]] +
  scale_fill_manual(values = as.vector(polychrome(15)[-c(1,2)]), na.value = "grey", name = NULL) +
  labs(x = "", y = "Order") +
  ylim(0,100)+
  theme(text = element_text(size = 22),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size =20))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/season_barplot_16s.png",width = 500,height = 500)
season.barplot.16s
dev.off() 

############### barplots 18S:

taxi <- tax.18s[["Order"]]
names(taxi) <- rownames(tax.18s)
relab.long <- melt(data.frame(use_relab.18s,otu = rownames(use_relab.18s)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["month"]][match(relab.long$variable,rownames(env))])
taxi <- taxi[match(relab.long$otu,names(taxi))]
agg.tax <- aggregate(relab.long$value,by = list(taxi),sum)
most.abu <- if(nrow(agg.tax)>12){
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)][1:12]
} else{
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)]
}
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)

colnames(relab.mean) <- c("otu","sub_cat","value")
use.tax <- tax.18s[match(relab.mean$otu,rownames(tax.18s)),]
relab.mean <- data.frame(relab.mean,use.tax)
relab.mean[["Order"]][!relab.mean[["Order"]]%in%c(most.abu,"unknown")] <- "Other"
relab.mean.agg  <- aggregate(relab.mean[,c("value")],by = list(relab.mean$sub_cat,relab.mean[["Order"]]),sum)
colnames(relab.mean.agg) <- c("Category","Taxonomy","Percent_cluster")
relab.mean.agg$Category <- factor(relab.mean.agg$Category, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
relab.mean.agg$Taxonomy <- gsub("_X","",relab.mean.agg$Taxonomy,fixed = T)
relab.mean.agg$Taxonomy <- gsub("Gammaproteobacteria_","",relab.mean.agg$Taxonomy,fixed = T)
relab.mean.agg$Taxonomy <- factor(relab.mean.agg$Taxonomy, levels = c(sort(unique(relab.mean.agg$Taxonomy[!(relab.mean.agg$Taxonomy%in%c("Other","Unknown"))])),"Other","Unknown"))

month.barplot.18s <- ggplot(relab.mean.agg, aes(x = Category, y = Percent_cluster, fill = Taxonomy)) +
  geom_bar(stat = "identity", color = "black") +
  themes[["Classic"]] +
  scale_fill_manual(values = as.vector(polychrome(15)[-c(1,2)]), na.value = "grey", name = NULL) +
  labs(x = "", y = "Order") +
  ylim(0,100)+
  theme(text = element_text(size = 22),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size =20))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/month_barplot_18s.png",width = 500,height = 500)
month.barplot.18s
dev.off() 

taxi <- tax.18s[["Order"]]
names(taxi) <- rownames(tax.18s)
relab.long <- melt(data.frame(use_relab.18s,otu = rownames(use_relab.18s)),id.vars = "otu")
relab.long <- relab.long %>% mutate(sub_cat = env[["Seasons"]][match(relab.long$variable,rownames(env))])
taxi <- taxi[match(relab.long$otu,names(taxi))]
agg.tax <- aggregate(relab.long$value,by = list(taxi),sum)
most.abu <- if(nrow(agg.tax)>12){
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)][1:12]
} else{
  agg.tax$Group.1[order(agg.tax$x,decreasing = T)]
}
relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)

colnames(relab.mean) <- c("otu","sub_cat","value")
use.tax <- tax.18s[match(relab.mean$otu,rownames(tax.18s)),]
relab.mean <- data.frame(relab.mean,use.tax)
relab.mean[["Order"]][!relab.mean[["Order"]]%in%c(most.abu,"unknown")] <- "Other"
relab.mean.agg  <- aggregate(relab.mean[,c("value")],by = list(relab.mean$sub_cat,relab.mean[["Order"]]),sum)
colnames(relab.mean.agg) <- c("Category","Taxonomy","Percent_cluster")
relab.mean.agg$Category <- factor(relab.mean.agg$Category, levels = c("Winter","Spring","Summer","Autumn"))
relab.mean.agg$Taxonomy <- factor(relab.mean.agg$Taxonomy, levels = c(sort(unique(relab.mean.agg$Taxonomy[!(relab.mean.agg$Taxonomy%in%c("Other","Unknown"))])),"Other","Unknown"))

season.barplot.18s <- ggplot(relab.mean.agg, aes(x = Category, y = Percent_cluster, fill = Taxonomy)) +
  geom_bar(stat = "identity", color = "black") +
  themes[["Classic"]] +
  scale_fill_manual(values = as.vector(polychrome(15)[-c(1,2)]), na.value = "grey", name = NULL) +
  labs(x = "", y = "Order") +
  ylim(0,100)+
  theme(text = element_text(size = 22),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size =20))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/season_barplot_18s.png",width = 500,height = 500)
season.barplot.18s
dev.off() 

############### Figure 2 of manuscript:

fig_2 <- plot_grid(month.barplot.16s,month.barplot.18s,plot.tax.king,plot.max.seasons,plot.full.temp,plot.full.chlo,align = "none",nrow = 2,
                       labels = c("A","B","C","D","E","F"),label_size = 35,label_colour = "black",label_x = c(-0.015,-0.015,0.01,0.01,0.01,0.01))

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/fig_2.png",width = 1600,height = 1000)
fig_2
dev.off()

############### Figure S4 of manuscript:

fig_s4 <-plot_grid(plot.full.daylength,plot.full.NH4,plot.full.NO2,plot.full.no3,plot.full.PO4,plot.full.si,plot.full.secchi,plot.full.sal,align = "none",nrow = 3,
                               labels = c("A","B","C","D","E","F","G","H"),label_size = 35,label_colour = "black")

png(filename = "~/Documents/spe_spe/Custom_figure_manuscript/Deutschmann/fig_s4.png",width = 1600,height = 1500)
fig_s4
dev.off()

