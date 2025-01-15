### Functions -------------------------------------------------------------


# Normalization function --------------------------------------------------
normal <- function(otu.table,type.norm,mean.thr = 0, max.thr = 0, occ.thr = 0,env.mat,sub.var,sub.cat,spars) {
  ### Get the right samples
  if(sub.var!= "Full network"&sub.cat!=""){
    if(sub.cat%in%unique(get(sub.var,env.mat))){
      use_table <- otu.table[,which(get(sub.var,env.mat)==sub.cat)]
    } else{
      use_table <- otu.table
    }
  } else{
    use_table <- otu.table
  }

  ### Filter data
    use_relab <- sweep(use_table, 2, colSums(use_table), FUN="/") * 100
    mmo <- NULL
    mmo$mean <- apply(use_relab, 1, mean, na.rm = T)
    mmo$max <- apply(use_relab, 1, max, na.rm = T)
    mmo$occ <- apply(use_relab > 0, 1, sum, na.rm = T)
    if(type.norm =="CLR"){
      use_table <- use_table[mmo$mean >= mean.thr & mmo$max >= max.thr & mmo$occ >= occ.thr,]
      set.seed(333)
      if(spars == "Random"){
      use_table[use_table==0] <- runif(sum(use_table==0),0.1,1)
      } else if(spars == "Constant"){
        use_table[use_table==0] <- 1
      }
      row.name <- rownames(use_table)
      normal_otu <- apply(use_table,2,function(x) log(x/mean(x)))
      rownames(normal_otu) <- row.name

    } else if(type.norm == "TSS"){
      normal_otu <- use_relab[mmo$mean >= mean.thr & mmo$max >= max.thr & mmo$occ >= occ.thr,]
      set.seed(333)
      if(spars == "Random"){
        normal_otu[normal_otu==0] <- runif(sum(normal_otu==0),0.1*min(normal_otu),min(normal_otu))
      }
      normal_otu <- sweep(normal_otu, 2, colSums(normal_otu), FUN="/") * 100
    } else if(type.norm == "SparCC"){
    normal_otu <- use_table[mmo$mean >= mean.thr & mmo$max >= max.thr & mmo$occ >= occ.thr,]
    set.seed(333)
    if(spars == "Random"){
      normal_otu[normal_otu==0] <- runif(sum(normal_otu==0),0.1,1)
    }
  }
  return(normal_otu)
}

select.tax <- function(tax.tab, mean.thr = 0, max.thr = 0, occ.thr = 0, mmo = mmo,relab.tab,env.mat,sub.var,sub.cat) {
  if(sub.var != "Full network"&sub.cat!=""){
    if(sub.cat%in%unique(get(sub.var,env.mat))){
      use_relab <- relab.tab[,which(get(sub.var,env.mat)==sub.cat)]
    } else{
      use_relab <- relab.tab
    }
  } else{
    use_relab <- relab.tab
  }
  use_relab <- sweep(use_relab, 2, colSums(use_relab), FUN="/") * 100
  mmo <- NULL
  mmo$mean <- apply(use_relab, 1, mean, na.rm = T)
  mmo$max <- apply(use_relab, 1, max, na.rm = T)
  mmo$occ <- apply(use_relab > 0, 1, sum, na.rm = T)
  use_tax <- tax.tab[mmo$mean >= mean.thr & mmo$max >= max.thr & mmo$occ >= occ.thr,]
  return(use_tax)
}

get.relabu <- function(otu.table, mean.thr = 0, max.thr = 0, occ.thr = 0,env.mat,sub.var,sub.cat) {
  if(sub.var != "Full network"&sub.cat!=""){
    if(sub.cat%in%unique(get(sub.var,env.mat))){
      use_relab <- otu.table[,which(get(sub.var,env.mat)==sub.cat)]
    } else{
      use_relab <- otu.table
    }
  } else{
    use_relab <- otu.table
  }
  use_relab <- sweep(use_relab, 2, colSums(use_relab), FUN="/") * 100
  mmo <- NULL
  mmo$mean <- apply(use_relab, 1, mean, na.rm = T)
  mmo$max <- apply(use_relab, 1, max, na.rm = T)
  mmo$occ <- apply(use_relab > 0, 1, sum, na.rm = T)
  use_otu <- use_relab[mmo$mean >= mean.thr & mmo$max >= max.thr & mmo$occ >= occ.thr,]
  set.seed(333)
  use_otu[use_otu==0] <- runif(sum(use_otu==0),0.1*min(use_otu),min(use_otu))
  use_otu <- sweep(use_otu, 2, colSums(use_otu), FUN="/") * 100
  return(use_otu)
}

select.rm <- function(relabu){
  rm.otu <- sqrt(rowMeans(relabu)) %>%
    rescale(to = c(1,12))
  return(rm.otu)
}

select.env <- function(env.tab,sub.var,sub.cat) {
  if(sub.var != "Full network"&sub.cat!=""){
    if(sub.cat%in%unique(get(sub.var,env.tab))){
      index <- get(sub.var,env.tab)==sub.cat
      index[is.na(index)]<- F
      use_env <- env.tab[index,]
    } else{
      use_env <- env.tab
    }
  } else{
    use_env <- env.tab
  }
  return(use_env)
}

env.plot.type <- function(env.var,env){
  if(env.var!=""&!is.null(env.var)){
    env.var <- get(env.var,env)
    if(is.numeric(env.var)){
      env.type <- "numeric"
    }
    if(is.character(env.var)|is.factor(env.var)){
      env.type <- "max"
    }
    if(length(unique(env.var))==1){
      env.type <- NULL
    }
  } else{
    env.type <- NULL
  }
  return(env.type)
}

# Correlation matrix based on abundance data ------------------------------
create.cor <- function(relab, type) {
  if(type%in%c("Spearman","Kendall","Pearson")){
    cor_mat <- cor(base::t(relab), method = tolower(type))
  } else if(type == "SparCC"){
    cor_mat <- t(sparcc(t(relab)))[[2]]
    rownames(cor_mat) <- rownames(relab)
    colnames(cor_mat) <- rownames(relab)
  }
  
  diag(cor_mat) <- NA
  return(cor_mat)
}

# Network construction ----------------------------------------------------
create.net <- function(matrix, cor.thr = 0.5) {
  otu_net <- graph.adjacency(matrix, weighted = TRUE, mode = "upper")
  otu_net <- delete.edges(otu_net, which(E(otu_net)$weight < cor.thr | E(otu_net)$weight %in% NA))
  otu_net <- igraph::simplify(otu_net, remove.multiple = T, remove.loops = T)
  return(otu_net)
}

# Calculate kmean clusters ------------------------------------------------
get.km <- function(method, matrix, no.of.clust = 3, nstart = 25, graph, use.seed = 1){
  set.seed(use.seed)
  if(method == "Cluster - kmeans"){
    diag(matrix) <- 1
    km.clust <- kmeans(matrix, no.of.clust, nstart)
    return(km.clust$cluster)
  }
}

# Calculate environmental parameter correlation ---------------------------
get.env.cor <- function(meta, use.env, matrix, type = "pairwise", plot.selection,net.norm){
  req(nrow(meta)==ncol(matrix)&use.env!="")
  if (plot.selection == "Environmental parameter"&use.env!="Not applicable"&use.env!="No numeric variables"){
    if(net.norm == "SparCC"){
      row.name <- rownames(matrix)
      normal_otu <- apply(matrix,2,function(x) log(x/mean(x)))
      rownames(matrix) <- row.name
    }
      if(is.numeric(meta[[use.env]])){
        env.cor <- as.numeric(cor(meta[[use.env]],base::t(matrix), use = type,method = "pearson"))
      }
      if(is.character(meta[[use.env]])|is.factor(meta[[use.env]])){
        relab.long <- melt(data.frame(matrix,otu = rownames(matrix)),id.vars = "otu")
        head(relab.long)
        relab.long <- relab.long %>% mutate(sub_cat = meta[[use.env]][match(relab.long$variable,rownames(meta))])
        relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)
        colnames(relab.mean) <- c("otu","sub_cat","value")
        max.cat <- sapply(unique(relab.mean$otu),function(x) relab.mean[relab.mean$otu==x,"sub_cat"][which.max(relab.mean[relab.mean$otu==x,"value"])])
        env.cor <- max.cat[match(rownames(matrix),names(max.cat))]
      }
      return(env.cor)
    }else{
      return(NULL)
    }
}

# Adjust taxonomy for plotting --------------------------------------------
get.tax <- function(tax, rank, graph, matrix, plot.selection = "tax") {
  if (plot.selection == "Taxonomy") {
    # Select taxonomic rank
    tax.sel <- tax %>% select(rank)
    tax.sel <- tibble::rownames_to_column(tax.sel, "name")
    
    # Merge graph nodes, mean abundances and taxonomy
    nodes <- data.frame(name = V(graph)$name)
    tax.flt <- tax.sel %>% filter(name %in% nodes$name)
    means <- data.frame(name = rownames(matrix),
                        meanAbundance = rowMeans(matrix))
    tax.mrg <- merge(means, tax.flt, by = "name", sort = FALSE)
    
    # Create list of top 11 most abundant taxonomies
    tax.top <- tax.mrg %>%
      group_by(.data[[rank]]) %>%
      summarise(sum = sum(meanAbundance)) %>%
      slice_max(sum, n = 15) %>%
      pull(rank)
    
    # Change taxonomy to "Other" for entries outside of the top 11
    final.tax <- case_when(tax.mrg[[rank]] %in% tax.top ~ tax.mrg[[rank]], TRUE ~ "Other")
    return(final.tax)
  }
}

# Convert network to tidy graph object ------------------------------------
tidy.net <- function(net, matrix, clust, env, tax,rm_otu) {
  req(!(nrow(matrix)!=length(env))|is.null(env))
  if(!is.null(tax)){
    tax <- factor(tax, levels = c(sort(unique(tax[!(tax%in%c("other","Other","unknown","Unknown"))])),"other","Other","unknown","Unknown"))
  }
  tidy.net <- net %>%
    as_tbl_graph()
  # Node size is proportional to the square root of the abundance of the corresponding OTU
  meanSize <- rm_otu 
  if(!is.null(env)){
    tidy.net %<>%
      activate(nodes) %>%
      mutate(nodeSize = meanSize,
             cluster = clust,
             envcor = env,
             taxo = tax,
             counts = rowMeans(matrix))
  } else{
    tidy.net %<>%
      activate(nodes) %>%
      mutate(nodeSize = meanSize,
             cluster = clust,
             taxo = tax,
             counts = rowMeans(matrix))
  }
  return(tidy.net)
}


# Plot themes -------------------------------------------------------------
themes <- list("Dark" = dark_theme_grey(),
               "Classic" = theme_bw(),
               "Minimal" = theme_minimal())

# Network plot function ---------------------------------------------------
plot.gg <- function(tidy.net, plot.selection, edge.width, edge.alpha, edge.strength,
                    seed, parameter, rank, theme, line.colour,isolated,clus,taxa,env.type,env,cor.env) {
  
  tidy.net
  
  lineCol <- list("Dark" = "white",
                  "Classic" = "black",
                  "Minimal" = "black")
  
  use_theme <- theme(axis.ticks = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     legend.title = element_text(size = 30),
                     legend.text = element_text(size = 30),
                     legend.position = "bottom",
                     legend.background = element_rect(fill='transparent'))#,
                     # legend.key.height = unit(0.25, 'cm'),
                     # legend.key.width = unit(0.3, "cm"))
  
  taxonomy <- apply(taxa,1,paste,collapse = ";")
  taxonomy <- gsub("'","",taxonomy)
  
  if(isolated=="No"){
    rem_nodes <- which(degree(tidy.net)==0)
    tidy.net <- delete_vertices(tidy.net,rem_nodes)
    if(length(rem_nodes)!=0){
      taxonomy <- taxonomy[-c(rem_nodes)]
    }
  }
  
  if (plot.selection == "Cluster - kmeans") {
    set.seed(seed)
    plot.clust <- ggraph(tidy.net, layout = 'fr') +
      geom_edge_arc0(
        width = edge.width,
        alpha = edge.alpha,
        strength = edge.strength,
        color = lineCol[[line.colour]]) +
      geom_point_interactive(aes(x,y,size = nodeSize,fill = as.factor(cluster),tooltip = taxonomy,data_id = taxonomy),pch=21,color = "black") +
      scale_size_continuous(range = c(1, 12), guide = 'none') +
      guides(fill = guide_legend("Cluster", override.aes = list(size = 15))) +
      scale_fill_manual(values = as.vector(polychrome(clus+3)[-c(1,2,3)])) +
      themes[[theme]] +
      use_theme
    plot.gi <- girafe(ggobj = plot.clust,width_svg = 15,height_svg = 15,
                      options= list(opts_hover_inv(css = "opacity:0.6;"),opts_hover(css = ''),opts_sizing(rescale = TRUE),
                                    opts_toolbar(hidden = c("lasso_select","lasso_deselect"))))
    return(plot.gi)
  }
  
  if (plot.selection == "Environmental parameter"&env.type =="numeric") {
    set.seed(seed)
    plot.env <- ggraph(tidy.net, layout = 'fr') +
      geom_edge_arc0(
        width = edge.width,
        alpha = edge.alpha,
        strength = edge.strength,
        color = lineCol[[line.colour]]) +
      geom_point_interactive(aes(x,y,size = nodeSize,fill = envcor,tooltip = taxonomy,data_id = taxonomy),pch=21,color = "black") +
      scale_size_continuous(range = c(1, 12), guide = 'none') +
      scale_fill_gradientn(colours = c("blue4", "dodgerblue", "#FFFFFF", "tomato", "red4"), limits = c(-1,1)) +
      themes[[theme]] +
      labs(fill = paste(parameter, "\ncorrelation")) +
      guides(fill = guide_colorbar(barheight = 5, barwidth = 15)) +
      use_theme+
      theme(legend.text = element_text(size=25))
    plot.gi <- girafe(ggobj = plot.env,width_svg = 15,height_svg = 15
                      ,options= list(opts_hover_inv(css = "opacity:0.6;"),opts_hover(css = ''),
                                     opts_sizing(rescale = TRUE),opts_toolbar(hidden = c("lasso_select","lasso_deselect"))))            
    return(plot.gi)
  }
  
  if (plot.selection == "Environmental parameter"&env.type == "max") {
    set.seed(seed)
      plot.env <- ggraph(tidy.net, layout = 'fr') +
        geom_edge_arc0(
          width = edge.width,
          alpha = edge.alpha,
          strength = edge.strength,
          color = lineCol[[line.colour]]) +
        geom_point_interactive(aes(x,y,size = nodeSize,fill = envcor,tooltip = taxonomy,data_id = taxonomy),pch=21,color = "black") +
        scale_size_continuous(range = c(1, 12), guide = 'none') +
        scale_fill_brewer(palette = "Paired", na.value = "grey") +
        themes[[theme]] +
        guides(fill = guide_legend(title = "Maximum mean\nrelative abundance in:",override.aes = list(size = 15),nrow=ceiling(length(unique(cor.env))/2),title.position = "top")) +
        use_theme+
        theme(legend.title.align = 0.5)
    plot.gi <- girafe(ggobj = plot.env,width_svg = 15,height_svg = 17,
                      options= list(opts_hover_inv(css = "opacity:0.6;"),opts_hover(css = ''),opts_sizing(rescale = TRUE),opts_toolbar(hidden = c("lasso_select","lasso_deselect"))))            
    return(plot.gi)
  }
  
  if (plot.selection == "Taxonomy") {
    set.seed(seed)
    leg.length <- min(ceiling(length(unique(vertex_attr(tidy.net, "taxo")))/2),8)
    
    plot.tax <- ggraph(tidy.net, layout = 'fr') +
      geom_edge_arc0(
        width = edge.width,
        alpha = edge.alpha,
        strength = edge.strength,
        color = lineCol[[line.colour]]) +
      geom_point_interactive(aes(x,y,size = nodeSize,fill = taxo,tooltip = taxonomy,data_id = taxonomy),pch=21,color = "black") +
      scale_size_continuous(range = c(1, 12), guide = 'none') +
      guides(fill = guide_legend(title = NULL, override.aes = list(size = 15),nrow = leg.length)) +
      scale_fill_manual(values = as.vector(polychrome(18)[-c(1,2)]), na.value = "grey") +
      themes[[theme]] +
      use_theme
    plot.gi <- girafe(ggobj = plot.tax,width_svg = 15,height_svg = 21,
                      options= list(opts_hover_inv(css = "opacity:0.6;"),opts_hover(css = ''),opts_sizing(rescale = TRUE),opts_toolbar(hidden = c("lasso_select","lasso_deselect"))))
    return(plot.gi)
  }

}

# Create taxonomy barplot -------------------------------------------------
tax.plot <- function(graph, option, yAxisTitle, legendTitle, yAxisData,env,env.var,env.type,plot.selection,env.main,use.relab,use.tax) {
# Scale mean abundances to 100% per cluster
  if(plot.selection == "Cluster - kmeans"){
    net_data <- graph %>%
      activate(nodes) %>%
      as_tibble() %>%
      group_by(cluster, taxo) %>%
      tally(counts) %>%
      mutate(freq = n / sum(n) * 100) %>%
      ungroup()
    net_data$cluster <- factor(net_data$cluster)
    colnames(net_data)[which(colnames(net_data)=="n")]<-"Percent_total"
    colnames(net_data)[which(colnames(net_data)=="freq")]<-"Percent_cluster"

    if(yAxisData=="Abundance over all clusters (%)"){
      taxplot <- ggplotly(
        ggplot(net_data, aes(x = cluster, y = Percent_total, fill = taxo)) +
          geom_bar(stat = "identity", color = "black") +
          themes[[option]] + 
          scale_fill_manual(values = as.vector(polychrome(18)[-c(1,2)]), na.value = "grey", name = NULL) +
          labs(x = "Cluster", y = legendTitle) +
          theme(text = element_text(size = 18))
      )
    }
    if(yAxisData%in%c("Abundance per cluster (%)","Abundance per category (%)")){
      taxplot <- ggplotly(
        ggplot(net_data, aes(x = cluster, y = Percent_cluster, fill = taxo)) +
          geom_bar(stat = "identity", color = "black") +
          themes[[option]] +
          scale_fill_manual(values = as.vector(polychrome(18)[-c(1,2)]), na.value = "grey",name = NULL) +
          labs(x = "Cluster", y = legendTitle) +
          theme(text = element_text(size = 18)), height = 700, width = 1250
      )
    }
  }
  if(plot.selection == "Taxonomy"|(plot.selection == "Environmental parameter"&env.type =="numeric")){
      net_data <- graph %>%
        activate(nodes) %>%
        as_tibble() %>%
        group_by(taxo) %>%
        tally(counts)
      
      net_data$x_var <- rep("All samples",nrow(net_data))
      colnames(net_data)[which(colnames(net_data)=="n")]<-"Percent_total"

      if(yAxisData=="Abundance over all clusters (%)"){
        taxplot <- ggplotly(
          ggplot(net_data, aes(x = x_var,y = Percent_total, fill = taxo)) +
            geom_bar(stat = "identity", color = "black") +
            themes[[option]] +
            scale_fill_manual(values = as.vector(polychrome(18)[-c(1,2)]), na.value = "grey", name = NULL) +
            labs(x = "All samples",y = legendTitle, fill = NULL) +
            theme(text = element_text(size = 18),axis.text.x = element_blank()), height = 700, width = 1250
        )
    }
      if(yAxisData%in%c("Abundance per cluster (%)","Abundance per category (%)")){
        taxplot <- ggplotly(
          ggplot(net_data, aes(x = x_var, y = Percent_total, fill = taxo)) +
            geom_bar(stat = "identity", color = "black") +
            themes[[option]] +
            scale_fill_manual(values = as.vector(polychrome(18)[-c(1,2)]), na.value = "grey", name = NULL) +
            labs(x= "All samples",y = legendTitle, fill = NULL) +
            theme(text = element_text(size = 18),axis.text.x = element_blank()), height = 700, width = 1250
        )
    }
  }
  
  if(plot.selection == "Environmental parameter"& env.type == "max"){
    taxi <- use.tax[[legendTitle]]
    names(taxi) <- rownames(use.tax)
    relab.long <- melt(data.frame(use.relab,otu = rownames(use.relab)),id.vars = "otu")
    relab.long <- relab.long %>% mutate(sub_cat = env[[env.var]][match(relab.long$variable,rownames(env))])
    taxi <- taxi[match(relab.long$otu,names(taxi))]
    agg.tax <- aggregate(relab.long$value,by = list(taxi),sum)
    agg.tax <- agg.tax[agg.tax$Group.1!="unknown",]
    most.abu <- if(nrow(agg.tax)>14){
      agg.tax$Group.1[order(agg.tax$x,decreasing = T)][1:14]
    } else{
      agg.tax$Group.1[order(agg.tax$x,decreasing = T)]
    }
    relab.mean <- aggregate(relab.long$value,by = list(relab.long$otu,relab.long$sub_cat),mean)

    colnames(relab.mean) <- c("otu","sub_cat","value")
    use_tax <- use.tax[match(relab.mean$otu,rownames(use.tax)),]
    relab.mean <- data.frame(relab.mean,use_tax)
    relab.mean[[legendTitle]][!relab.mean[[legendTitle]]%in%c(most.abu,"unknown")] <- "Other"
    relab.mean.agg  <- aggregate(relab.mean[,c("value")],by = list(relab.mean$sub_cat,relab.mean[[legendTitle]]),sum)
    colnames(relab.mean.agg) <- c("Category","Taxonomy","Percent_cluster")
    
      taxplot <- ggplotly(
        ggplot(relab.mean.agg, aes(x = Category, y = Percent_cluster, fill = Taxonomy)) +
          geom_bar(stat = "identity", color = "black") +
          themes[[option]] +
          scale_fill_manual(values = as.vector(polychrome(18)[-c(1,2)]), na.value = "grey", name = NULL) +
          labs(x = "", y = legendTitle) +
          theme(text = element_text(size = 18),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size =20)), height = 900, width = 1250
      ) 
  }
  return(taxplot)
}

# Create corplot ----------------------------------------------------------
cor.plot <- function(graph, parameter, option, Correlation,clus,env,env.var,env.type,plot.selection,env.main) {
  if(!is.null(Correlation)){
    if(plot.selection == "Cluster - kmeans"){
      cor.frame <- graph %>%
        as_tibble()
      cor.frame$cluster <- factor(cor.frame$cluster)
      colnames(cor.frame)[colnames(cor.frame)=="taxo"] <- "Taxonomy"
      cor.plot<-ggplotly(
          ggplot(cor.frame,aes(x = cluster, y = Correlation)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(col = Taxonomy),size = 0.7) +
          scale_color_manual(values = as.vector(polychrome(18)[-c(1,2)])) +
          geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth = .5) +
          themes[[option]] +
          labs(x = "Cluster", y = paste(parameter, "\ncorrelation"), fill = "Cluster") +
          theme(text=element_text(size=18)), height = 600, width = 1250)

      
      for(i in 1:length(cor.plot$x$data)){
        if("outliercolor"%in%names(cor.plot$x$data[[i]]$marker)){
          cor.plot$x$data[[i]]$marker$opacity = 0
        }
      }
      
    }
    if(plot.selection == "Taxonomy"|plot.selection == "Environmental parameter"&env.type =="numeric"){
      cor.frame <- graph %>%
        as_tibble()
      colnames(cor.frame)[colnames(cor.frame)=="taxo"] <- "Taxonomy"
      cor.plot<-ggplotly(
          ggplot(cor.frame, aes(x="",y = Correlation)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(col = Taxonomy),size = 0.7) +
          scale_color_manual(values = as.vector(polychrome(18)[-c(1,2)])) +
          geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth = .5) +
          themes[[option]] +
          labs(x = "All clades", y = paste(parameter, "\ncorrelation")) +
          theme(text=element_text(size=18),axis.text.x = element_blank(),axis.ticks.x = element_blank()), height = 600, width = 1250)

      for(i in 1:length(cor.plot$x$data)){
        if("outliercolor"%in%names(cor.plot$x$data[[i]]$marker)){
          cor.plot$x$data[[i]]$marker$opacity = 0
        }
      }
    }
    if(plot.selection == "Environmental parameter"& env.type == "max"){
      if(env.type == "max"){
        cols.nr = length(unique(env.main))
      }
      cor.frame <- graph %>%
        as_tibble()
      colnames(cor.frame)[colnames(cor.frame)=="taxo"] <- "Taxonomy"
      cor.plot<-ggplotly(
          ggplot(cor.frame,aes(x = as.factor(env.main), y = Correlation)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(col = Taxonomy),size = 0.7) +
          scale_color_manual(values = as.vector(polychrome(18)[-c(1,2)])) +
          geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth = .5) +
          themes[[option]] +
          labs(x = "Environment", y = paste(parameter, "\ncorrelation"), fill = "Environment") +
          theme(text=element_text(size=18),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 20)), height = 600, width = 1250)
      
      for(i in 1:length(cor.plot$x$data)){
        if("outliercolor"%in%names(cor.plot$x$data[[i]]$marker)){
          cor.plot$x$data[[i]]$marker$opacity = 0
        }
      }
    }
    return(cor.plot)
    } else{
    return(NULL)
    }
}


