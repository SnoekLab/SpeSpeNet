###################### This scripts generates figure S5 of the paper titled "SpeSpeNet: An interactive and user-friendly tool 
###################### to create and explore microbial correlation networks." The data used is from Hauptfeld et al 2022 
###################### https://doi.org/10.1016/j.watres.2022.118767

library(NetCoMi)
library(phyloseq)
library(pals)
library(tidyr)
library(plotrix)

load("NW_data/Griftpark/grift.phy.out")

net_grift <- netConstruct(grift.phy,taxRank = "Genus",dataType = "counts",filtTax = c("numbSamp","relFreq"),
                          filtTaxPar = list(numbSamp = 5,relFreq = 0.0005),filtSamp = "none",measure = "sparcc",
                          zeroMethod = "pseudoZO",zeroPar = list(pseudoZO = 0.1),sparsMethod = "threshold",thresh = 0.5,
                          seed = 123)

net_grift_ana <- netAnalyze(net_grift)

set.seed(4)

png(filename = "~/Documents/spe_spe/figures_manuscript/NetCoMi/clusters.png",width = 700,height = 500)
plot(net_grift_ana,layout = "layout_with_fr",nodeColor = "cluster",nodeSize = "TSS",negCol = "white",nodeSizeSpread =6,
     featVecCol = orders,labels = F,edgeWidth = 0.1,posCol = "grey",highlightHubs = F,
     nodeTransp = 10)
legend(0.3, 1, cex = 1.5, pt.cex = 2, title = "Cluster", 
       legend=factor(c(1,2)), col = c("cyan","red"), bty = "n", pch = 16) 
dev.off()

taxtab <- as(tax_table(grift.phy), "matrix")

orders <- taxtab[, "Order"]
names(orders) <- taxtab[, "Order"]

orders <- as.factor(ifelse(orders%in%c("Burkholderiales","Campylobacterales","Deltaproteobacteria","Desulfobacterales",
                                       "Pseudomonadales","Rhodocyclales","Sphingomonadales","Xanthomonadales","Unknown"),
                           orders,"Other"))
names(orders) <- taxtab[,6]

ordercol <- as.vector(c(polychrome(12)[-c(1,2)],"grey"))

set.seed(4)

png(filename = "~/Documents/spe_spe/figures_manuscript/NetCoMi/orders.png",width = 700,height = 500)
plot(net_grift_ana,layout = "layout_with_fr",nodeColor = "feature",nodeSize = "TSS",negCol = "white",nodeSizeSpread =6,
     featVecCol = orders,colorVec =  ordercol,labels = F,edgeWidth = 0.1,posCol = "grey",highlightHubs = F,
     nodeTransp = 10)
legend(0.3, 1, cex = 1.2, pt.cex = 2, title = "Order", 
       legend=levels(orders), col = ordercol, bty = "n", pch = 16) 
dev.off()


otu <- data.frame(as(otu_table(grift.phy), "matrix"))
env <- as(sample_data(grift.phy), "matrix")
env <- data.frame(env)

otu[otu==0] <- runif(sum(otu==0),0.1,1)
row.name <- rownames(otu)
otu <- apply(otu,2,function(x)log(x/mean(x)))
rownames(otu) <- row.name
otu <-data.frame(otu)

cor.vec <- as.numeric(cor(as.numeric(env[["O2.levels"]]),base::t(otu), use = "pairwise",method = "pearson"))

names(cor.vec) <- taxtab[,6]

ramp <- colorRamp(c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"))
use.col <- rgb( ramp(seq(0, 1, length.out = 256)), maxColorValue = 256)

brks <- seq(min(cor.vec), max(cor.vec), length.out = 256)
grp <- cut(cor.vec, breaks = brks, include.lowest = TRUE)
col.vec <- use.col[grp]
names(col.vec) <- names(cor.vec)

set.seed(4)
png(filename = "~/Documents/spe_spe/figures_manuscript/NetCoMi/O2_cor.png",width = 700,height = 500)
plot(net_grift_ana,layout = "layout_with_fr",nodeColor = "colorVec",nodeSize = "TSS",negCol = "white",nodeSizeSpread =6,
     featVecCol = cor.vec,colorVec =  col.vec,labels = F,edgeWidth = 0.1,posCol = "grey",highlightHubs = F,
     nodeTransp = 10)
color.legend(0.65,-0.3,0.3,0.6,c(-0.1,-0.05,0,0.05,1),use.col,gradient="y",align = "rb")
dev.off()

