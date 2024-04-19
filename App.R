

source("nw_functions.R")

options(shiny.maxRequestSize = 40 * 1024^2)
options(getClass.msg=FALSE)

# Define UI ---------------------------------------------------------------

ui <- fluidPage(titlePanel(
  tagList(
    img(src = "SpeSpe-logos_white.png", height = 80),
    "v0.1.8",
    img(src = "uu_logo_transparent.png", height = 80,
        style = {"position:absolute;right:15px;z-index:1000000;"}
    )
  )
),
theme = shinytheme("cyborg"),

tabsetPanel(
  
  # Network tab -------------------------------------------------------------
  
  tabPanel("Network",
           sidebarLayout(
             sidebarPanel(width = 3,
                          selectInput("indata", "Select input type", choices = c(".txt file", "Phyloseq", "MGnify")),
                          conditionalPanel(
                            condition = "input.indata == '.txt file'",
                            selectInput("normalize.custom",
                                        label = tags$span(
                                          "Are OTU counts normalized?"),
                                        choices = c("No - raw read counts", "Yes - percentages (0 to 100%)", "Yes - fractions (0 to 1)")),
                                        selectInput("agglomerate.custom",label = "Aggregate OTUs/ASVs at genus level?", choices = c("Yes", "No")),
                            fileInput("otu.file", "Submit OTU table", multiple = TRUE,
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".csv",".txt"), buttonLabel = "Browse..."),
                            fileInput("tax.file", "Submit Tax", multiple = TRUE,
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".csv",".txt"), buttonLabel = "Browse..."),
                            fileInput("env.file", "Submit Env table", multiple = TRUE,
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".csv",".txt"), buttonLabel = "Browse..."),
                          ),
                          conditionalPanel(
                            condition = "input.indata == 'Phyloseq'",
                            selectInput("normalize.phylo",
                                        label = tags$span(
                                          "Are OTU counts normalized?"),
                                        choices = c("No - raw read counts", "Yes - percentages (0 to 100%)", "Yes - fractions (0 to 1)")
                                        ),
                            selectInput("agglomerate.phylo",label = "Aggregate OTUs/ASVs at genus level?", choices = c("Yes", "No")),
                            fileInput("phy.file", "Submit Phyloseq .rds file", buttonLabel = "Browse...")
                          ),
                          conditionalPanel(
                            condition = "input.indata == 'MGnify'",
                            textInput("mgnify.acc",label = tags$span("Submit MGnify study accession",bsButton("mgInfo", label = "",icon = icon("info"),
                                                                                                               style = "info", size = "extra-small")), 
                                      placeholder = "MGYS00001447", value = "MGYS00001447"),
                            selectInput("agglomerate.mgnify",label = "Aggregate OTUs/ASVs at genus level?", choices = c("Yes", "No")),
                            numericInput("max.sample", "Number of samples to download", value = 20, min = 1, max = 200, step = 1),
                            selectInput("functions","include functional data?",choices = c("Yes","No"))
                          ),
                          actionButton("load.btn", "Load data", icon = icon("fas fa-send", verify_fa = FALSE)),
                          tags$hr(),
                          h4("Network construction settings"),
                          sliderInput('cor.thr', label = "Correlation threshold", 0.5, min = 0.1, max = 0.99, step = 0.01),
                          numericInput("occ", label = "Occurrence threshold", value = 5, step = 1),
                          numericInput("max.ab", label = "Minimal abundance threshold (%)", value = 0.1, step = 0.05, max = 100),
                          h4("Aesthethic mapping"),
                          fluidRow(
                            column(width = 8,
                                   selectInput("colors", "Color by", choices = c("Cluster - kmeans","Environmental parameter", "Taxonomy"))),
                            column(width = 4,
                                   numericInput("set.seed", "Set seed", value = 1, min = 1, step = 1))
                          ),
                          fluidRow(
                            column(width = 8,
                                   conditionalPanel(
                                     condition = "input.colors == 'Cluster - kmeans'",
                                     numericInput("nclust", "n clusters", 3, min = 1, max = 20, step = 1)
                                   ),
                                   conditionalPanel(
                                     condition = "input.colors == 'Environmental parameter'",
                                     selectInput("env.network", "Select environmental parameter", choices = NULL)
                                   ),
                                   conditionalPanel(
                                     condition = "input.colors == 'Taxonomy'",
                                     selectInput("tax.network", "Select taxonomic rank", choices = NULL)
                                   )
                            ),
                            fluidRow(
                              column(width = 6,
                                     selectInput("sub.variable","Select variable for sub network", choices = c("Full network"),selected = "Full network")
                                     ),
                              column(width = 6,
                                     selectInput("sub.category","Select category for sub network", choices = NULL)
                                     )
                            ),
                            column(width = 4,
                                   selectInput("net.theme", "Plot theme", choices = names(themes))
                            )
                          ),
                          fluidRow(
                            column(width = 3,
                                   numericInput("edge.arc", "Arc strength", value = 0, min = -1, max = 1, step = 0.1)),
                            column(width = 3,
                                   numericInput("edge.width", "Edge Width", value = 0.1, min = 0.1, max = 2, step = 0.1)),
                            column(width = 3,
                                   numericInput("edge.alpha", "Edge Alpha", value = 1, min = 0.1, max = 1, step = 0.1)),
                            column(width = 3,
                                   selectInput("iso.network", "Plot isolated nodes", choices = c("Yes", "No"))),
                            column(width = 4,
                                   selectInput("cor.method","Correlation method", choices = c("spearman","pearson","kendall")))
                          )),
             mainPanel(girafeOutput("network", height = 1100, width = 1100) %>% withSpinner(type = 6),
                       htmlOutput("netstats"))
           )),
  
  # Summary tab -------------------------------------------------------------
  
  tabPanel("Summary",
           sidebarLayout(
             sidebarPanel(width = 3,
                          selectInput("tax.summary", "Select taxonomic rank", choices = NULL),
                          selectInput("env.summary", "Select environmental parameter", choices = NULL),
                          selectInput("tax.axis", "Barplot y-axis", choices =c("Abundance per cluster (%)", "Abundance over all clusters (%)")),
                          selectInput("sum.theme", "Select plot theme", choices = names(themes))),
             mainPanel(plotlyOutput("sumplot") %>% withSpinner(type = 6),height = "2000px"))),
  
  # Raw data downloads ------------------------------------------------------
  navbarMenu(
    "Raw data",
    tabPanel("Network data",
             sidebarLayout(
               sidebarPanel(width = 3,
                            h4("Download raw network data"),
                            selectInput("preview", "Select preview", choices = c("Node data", "Edge data", "Similarity matrix")),
                            downloadButton("dl.matrix", "Similarity matrix"),
                            downloadButton("dl.network", "Edge data"),
                            downloadButton("dl.aesthetic", "Node mappings")),
               mainPanel(dataTableOutput("rawtables"))
             ))
  )
)
)

server <- function(input, output, session) {
    
  options(shiny.maxRequestSize=300*1024^2)
  
  # Load data ---------------------------------------------------------------
  
  # Convert MGnify analyses to phyloseq object
  
  mgnify_phylo <- eventReactive(input$load.btn, {
    if (input$indata == "MGnify"){
      mgclnt <- MgnifyClient(usecache = T, cache_dir = '/tmp/MGnify_cache')
      accession <- input$mgnify.acc
      # Gather analyses accessions
      mgnify_analyses <- searchAnalysis(mgclnt, accession,type = "studies")  
      # Adjust max analyses able to be pulled
      if(as.numeric(input$max.sample) > length(mgnify_analyses)){
        mgn_phylo <- getResult(mgclnt, mgnify_analyses, get.func = F,get.taxa = T,output = "phyloseq"
                               ,bulk_dl = T,get.tree = F,use.cache = T)
      } 
      if(as.numeric(input$max.sample)<= length(mgnify_analyses)){
        mgn_phylo <- getResult(mgclnt,mgnify_analyses[1:as.numeric(input$max.sample)], get.func = F, get.taxa = T,output = "phyloseq"
                                            ,bulk_dl = T,get.tree = F,use.cache = T)
      }
      # Store sample, tax and otu data in separate objects
      sam_phylo <- mgn_phylo@sam_data %>% sample_data()%>% data.frame()
      tax_phylo <- mgn_phylo@tax_table %>% tax_table()%>% data.frame()
      otu_phylo <- mgn_phylo@otu_table %>% otu_table() %>% data.frame()
      otu_phylo <- normal(data = otu_phylo, normalization = "No - raw read counts")
      return(list(sam = sam_phylo, tax = tax_phylo, otu = otu_phylo))
    }
  })
  
  ## Env file ---------------------------------------------------------------
  env_data <- eventReactive(input$load.btn, {
    if(input$indata == ".txt file"){
      if(is.null(input$env.file)){
        return(NULL)
      }
      else{
        env <- read.delim(input$env.file$datapath, header = T, sep = "\t",row.names = 1,check.names = F)
        env_num <- env %>%
          select(where(is.numeric))
        env_num <- env_num[,which(apply(env_num,2,function(x)length(unique(x))>1))]
        env_fac <- env %>% 
          select(where(is.factor)|where(is.character))
        env_fac <- env_fac[which(apply(env_fac,2,function(x) length(unique(x))>1))]
        env_sub <- env_fac[which(apply(env_fac,2,function(x) max(table(x),na.rm = T)>7))]
        updateSelectInput(session, "sub.variable", label = "Select variable for sub network", choices = c("Full network",colnames(env_sub)), selected = "Full network")
        updateSelectInput(session, "env.network", label = "Select environmental parameter", choices  = colnames(env_num), selected = colnames(env_num)[1])
        updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = colnames(env_num), selected = colnames(env_num)[1])
        return(env)
      }
    }
    if(input$indata == "Phyloseq"){
      phy_env <- readRDS(input$phy.file$datapath)
      phy_env <- phy_env %>%
        sample_data() %>%
        data.frame()
      phy_fac <- phy_env %>% 
        select(where(is.factor)|where(is.character))
      phy_fac <- phy_fac[which(apply(phy_fac,2,function(x)length(unique(x))>1))]
      phy_sub <- phy_fac[which(apply(phy_fac,2,function(x) max(table(x),na.rm = T)>7))]
      if(length(apply(phy_env,2,function(x)is.numeric(x)))>0){
        phy_num <- phy_env %>%
          select(where(is.numeric))
        phy_num <- phy_num[,which(apply(phy_num,2,function(x)length(unique(x))>1))]
        updateSelectInput(session, "sub.variable", label = "Select variable for sub network", choices = c("Full network",colnames(phy_sub)), selected = "Full network")
        updateSelectInput(session, "env.network", label = "Select environmental parameter", choices = colnames(phy_num), selected = colnames(phy_num)[1])
        updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices = colnames(phy_num), selected = colnames(phy_num)[1])
      }
      return(phy_env)
    }
    if(input$indata == "MGnify"){
      mgnify_sam <- mgnify_phylo()[['sam']]
      mgnify_sam <- type.convert(mgnify_sam, as.is = T)
      mgnify_num <- mgnify_sam %>%
        select(where(is.numeric))
      mgnify_num <- mgnify_num[,which(apply(mgnify_num,2,function(x)length(unique(x))>1))]
      mgnify_fac <- mgnify_sam %>% 
        select(where(is.factor)|where(is.character))
      mgnify_fac <- mgnify_fac[which(apply(mgnify_fac,2,function(x) length(unique(x))>1))]
      mgnify_sub <- mgnify_fac[which(apply(mgnify_fac,2,function(x) max(table(x),na.rm = T)>7))]
      updateSelectInput(session, "sub.variable", label = "Select variable for sub network", choices = c("Full network",colnames(mgnify_sub)), selected = "Full network")
      updateSelectInput(session, "env.network", label = "Select environmental parameter", choices  = colnames(mgnify_num), selected = colnames(mgnify_num)[1])
      updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = colnames(mgnify_num), selected = colnames(mgnify_num)[1])
      return(mgnify_sam)
    }
  })
  
  
  eventReactive(input$load.btn, {
    meta <- env_data()
    env.vars <- sapply(meta,function(x) (length(unique(x))>1&(length(unique(x))<13|is.numeric(x))))
    env.vars <- colnames(meta)[env.vars]
    if(input$indata == ".txt file") {
      updateSelectInput(session, "tax.network", label = "Select taxonomic rank", choices  = colnames(tax_data()$tax))
      updateSelectInput(session, "tax.summary", label = "Select taxonomic rank", choices  = colnames(tax_data()$tax))
      updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = colnames(env_data()$env_num))
      updateSelectInput(session, "env.network", label = "Select environmental parameter", choices  = env.vars)

    }
    if(input$indata == "Phyloseq") {
      updateSelectInput(session, "tax.network", label = "Select taxonomic rank", choices  = colnames(tax_data()$phy_tax))
      updateSelectInput(session, "tax.summary", label = "Select taxonomic rank", choices  = colnames(tax_data()$phy_tax))
      updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = colnames(env_data()$phy_num))
      updateSelectInput(session, "env.network", label = "Select environmental parameter", choices  = env.vars)
    }
    if(input$indata == "MGnify") {
      updateSelectInput(session, "tax.network", label = "Select taxonomic rank", choices  = colnames(tax_data()$mgnify_tax))
      updateSelectInput(session, "tax.summary", label = "Select taxonomic rank", choices  = colnames(tax_data()$mgnify_tax))
      updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = colnames(env_data()$mgnify_num))
      updateSelectInput(session, "env.network", label = "Select environmental parameter", choices  = env.vars)
    }
  })
  
  ## OTU file ---------------------------------------------------------------
  otu_data <- eventReactive(input$load.btn, {
    if(input$indata == ".txt file"){
      env <- read.delim(input$env.file$datapath, header = T, sep = "\t",row.names = 1,check.names = F)
      relab <- read.delim(input$otu.file$datapath, header = T, sep = "\t",row.names = 1,check.names = F)
      tax <- read.delim(input$tax.file$datapath, header = T, sep = "\t",row.names = 1,check.names = F)
      validate(
        need(isTRUE(all.equal(colnames(relab),rownames(env))), "column names OTU file don't match row names metadata file"),
        need(isTRUE(all.equal(rownames(relab),rownames(tax))), "Row names OTU file don't match row names taxonomy file"))
      relab <- normal(data = relab, normalization = input$normalize.custom)
      if(input$agglomerate.custom == "Yes"){
        tax <- read.delim(input$tax.file$datapath, header = T, sep = "\t",row.names = 1,check.names = F)
        tax[is.na(tax)] <- "Unknown"
        tax[tax == "NA"] <- "Unknown"
        colnames(tax) <- toTitleCase(colnames(tax))
        relab <- aggregate(relab,tax[,tolower(colnames(tax))!="species"],sum)
        relab <- relab %>%
          select(where(is.numeric))
      } 
      return(relab)
    }
    if(input$indata == "Phyloseq"){
      phylo <- readRDS(input$phy.file$datapath)
      phy_otu <- phylo %>%
        otu_table() %>%
        as.data.frame()
      nr.otu <- nrow(phylo@tax_table)
      if(ncol(phy_otu)==nr.otu){
        phy_otu <- t(phy_otu)
      }
      new_otu <- normal(data = phy_otu, normalization = input$normalize.phylo)
      if(input$agglomerate.phylo == "Yes"){
        tax <- phylo %>%
          tax_table() %>%
          as.data.frame()
        tax[is.na(tax)] <- "Unknown"
        tax[tax == "NA"] <- "Unknown"
        colnames(tax) <- toTitleCase(colnames(tax))
        new_otu <- aggregate(new_otu,tax[,colnames(tax)!="Species"],sum)
        new_otu <- new_otu %>%
          select(where(is.numeric))
      }
      return(new_otu)
    }
    if(input$indata == "MGnify"){
      mgnify_otu <- mgnify_phylo()[['otu']]
      if(input$agglomerate.mgnify == "Yes"){
        tax <- mgnify_phylo()[['tax']]
        tax[is.na(tax)] <- "Unknown"
        tax[tax == "NA"] <- "Unknown"
        colnames(tax) <- toTitleCase(colnames(tax))
        mgnify_otu <- aggregate(mgnify_otu,tax[,colnames(tax)!="Species"],sum)
        mgnify_otu <- mgnify_otu %>%
          select(where(is.numeric))
      }
      return(mgnify_otu)
    }
  })
  
  ## Tax file --------------------------------------------------------------
  tax_data <- eventReactive(input$load.btn, {
    if(input$indata == ".txt file") {
      if(is.null(input$tax.file)){
        return(NULL)
      }
      else{
        tax <- read.delim(input$tax.file$datapath, header = T, sep = "\t",row.names = 1,check.names = F)
        tax[is.na(tax)] <- "Unknown"
        tax[tax == "NA"] <- "Unknown"
        colnames(tax) <- toTitleCase(colnames(tax))
        if(input$agglomerate.custom == "Yes"){
          relab <- read.delim(input$otu.file$datapath, header = T, sep = "\t",row.names = 1,check.names = F)
          relab <- normal(data = relab, normalization = input$normalize.custom)
          relab <- aggregate(relab,tax[,colnames(tax)!="Species"],sum)
          tax <- relab %>% select(where(is.character))
        } 
        updateSelectInput(session, "tax.network", label = "Select taxonomic rank", choices  = colnames(tax), selected = colnames(tax)[2])
        updateSelectInput(session, "tax.summary", label = "Select taxonomic rank", choices  = colnames(tax), selected = colnames(tax)[2])
      }
      return(tax)
    }
    if(input$indata == "Phyloseq") {
      phylo <- readRDS(input$phy.file$datapath)
      phy_tax <- phylo %>%
        tax_table() %>%
        as.data.frame()
      phy_tax[is.na(phy_tax)] <- "Unknown"
      phy_tax[phy_tax == "NA"] <- "Unknown"
      colnames(phy_tax) <- toTitleCase(colnames(phy_tax))
      if(input$agglomerate.phylo == "Yes"){
        phy_otu <- phylo %>%
          otu_table() %>%
          as.data.frame()
        nr.otu <- nrow(phylo@tax_table)
        if(ncol(phy_otu)==nr.otu){
          phy_otu <- t(phy_otu)
        }
        new_otu <- normal(data = phy_otu, normalization = input$normalize.phylo)
        new_otu <- aggregate(new_otu,phy_tax[,colnames(phy_tax)!="Species"],sum)
        phy_tax <- new_otu %>% select(where(is.character))
      } 
      updateSelectInput(session, "tax.network", label = "Select taxonomic rank", choices  = colnames(phy_tax), selected = colnames(phy_tax)[2])
      updateSelectInput(session, "tax.summary", label = "Select taxonomic rank", choices  = colnames(phy_tax), selected = colnames(phy_tax)[2])
      return(phy_tax)
    }
    if(input$indata == "MGnify"){
      mgnify_tax <- mgnify_phylo()[['tax']]
      mgnify_tax[is.na(mgnify_tax)] <- "Unknown"
      mgnify_tax[mgnify_tax == "NA"] <- "Unknown"
      colnames(mgnify_tax) <- toTitleCase(colnames(mgnify_tax))
      if(input$agglomerate.mgnify == "Yes"){
        mgnify_otu <- mgnify_phylo()[['otu']]
        mgnify_otu <- aggregate(mgnify_otu,mgnify_tax[,colnames(mgnify_tax)!="Species"],sum)
        mgnify_tax <- mgnify_otu %>% select(where(is.character))
      } 
      updateSelectInput(session, "tax.network", label = "Select taxonomic rank", choices  = colnames(mgnify_tax), selected = colnames(mgnify_tax)[2])
      updateSelectInput(session, "tax.summary", label = "Select taxonomic rank", choices  = colnames(mgnify_tax), selected = colnames(mgnify_tax)[2])
      return(mgnify_tax)
    }
  })
  
  observeEvent(input$sub.variable, {
    if(input$sub.variable!="Full network"){
      cat.choice <- table(get(input$sub.variable,env_data()))
      cat.choice <- names(cat.choice)[cat.choice>7]
      updateSelectInput(session,"sub.category", label = "Select category for sub network", choices = cat.choice,selected = cat.choice[1])
    }
    if(input$sub.variable=="Full network"){
      updateSelectInput(session,"sub.category", label = "Select category for sub network", choices = c("Show full network"), selected = "Show full network")
    }
  })
  

    use.env <- reactive({
    sub.var <- isolate(input$sub.variable)
    select.env(env.tab = env_data(),sub.var=sub.var, sub.cat = input$sub.category)})
  
  observeEvent(use.env(),{
    meta <- use.env()
    sub.var <- input$sub.variable
    env.vars <- sapply(meta,function(x) (length(unique(x))>1&(length(unique(x))<13|is.numeric(x))))
    env.vars <- colnames(meta)[env.vars]
    env.vars <- env.vars[env.vars!=sub.var]
    if(!input$env.network%in%env.vars){
      updateSelectInput(session, "env.network", label = "Select environmental parameter", choices  = env.vars,selected = env.vars[1])
    } else{
      updateSelectInput(session, "env.network", label = "Select environmental parameter", choices  = env.vars,selected = input$env.network)
    }
      
  })
  
  use.otu <- reactive({
    sub.var <- isolate(input$sub.variable)
    select.genera(relab.tab = otu_data() , mean.thr = 0,
                                     max.thr = input$max.ab, occ.thr = input$occ,
                                     env.mat = env_data(), sub.var = sub.var, sub.cat = input$sub.category)})
  
  use.tax <- reactive({
    sub.var <- isolate(input$sub.variable)
    select.tax(tax.tab = tax_data() , mean.thr = 0,
                                     max.thr = input$max.ab, occ.thr = input$occ, relab.tab = otu_data(),
                                  env.mat=env_data(),sub.var = sub.var,sub.cat = input$sub.category)})

  
  env.col <- reactive({
    env.plot.type(env.var = input$env.network,env = isolate(use.env()))
  })
  
  observeEvent(list(env.col(),input$colors),{
    req(!is.null(env.col()))
    type.var <- input$colors
    type.env <- env.col()
    if(type.env=="max"&type.var=="Environmental parameter"){
      updateSelectInput(session, "tax.axis", label = "Barplot y-axis", choices =c("Abundance per category (%)"))
      updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = "Not applicable", selected = "Not applicable")

    } else{
      meta <- use.env()
      env_num <- meta %>%
        select(where(is.numeric))
      env_num <- env_num[,which(apply(env_num,2,function(x)length(unique(x))>1))]
      updateSelectInput(session, "tax.axis", "Barplot y-axis", choices =c("Abundance per cluster (%)", "Abundance over all clusters (%)"))
      if(!is.null(dim(env_num))){
        updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = colnames(env_num), selected = colnames(env_num)[1])
      } else{
        updateSelectInput(session, "env.summary", label = "Select environmental parameter", choices  = "No numeric variables", selected = "No numeric variables")
      }
    }
  })
  
  ## Correlation matrix ------------------------------------------------------
  cor.mat <- reactive({
    create.cor(relab = use.otu(), type = input$cor.method)})
  
  ## Calculate environmental correlation -------------------------------------
  env.cor <- reactive({
  req(input$env.network)
    input$occ
    input$max.ab
    get.env.cor(meta = use.env(), use.env = input$env.network, matrix = isolate(use.otu()), plot.selection = input$colors)})
  
  ## Graph construction ------------------------------------------------------
  net.obj <- reactive({
    create.net(matrix = cor.mat(), cor.thr = input$cor.thr)})
  
  ## Calculate kmeans clusters ------------------------------------------------
  clust.obj <- reactive({get.km(matrix = cor.mat(), no.of.clust = input$nclust, graph = net.obj(), method = input$colors)})
  
  ## Create taxonomy for plotting --------------------------------------------
  tax.obj <- reactive({get.tax(tax_data(), input$tax.network, net.obj(), use.otu(), plot.selection = input$colors)})
  
  ## Convert to tidy graph ---------------------------------------------------
  net.tbl <- reactive({
    tidy.net(net = net.obj(), matrix = isolate(use.otu()), clust = clust.obj(), env = env.cor(), tax = tax.obj())})
  
  ## Plot network ------------------------------------------------------------
  spe.plt <- reactive({
    plot.gg(
      tidy.net = net.tbl(),
      plot.selection = isolate(input$colors),
      seed = input$set.seed,
      parameter = isolate(input$env.network),
      rank = isolate(input$tax.network),
      edge.alpha = input$edge.alpha,
      edge.width = input$edge.width,
      edge.strength = input$edge.arc,
      theme = input$net.theme,
      line.colour = input$net.theme,
      isolated = input$iso.network,
      clus = isolate(input$nclust),
      taxa = isolate(use.tax()),
      env.type = isolate(env.col()),
      env = isolate(use.env()),
      cor.env = isolate(env.cor())
    )
  })
  
  # Plot network ------------------------------------------------------------
  output$network <- renderGirafe({
    validate(
      need(input$load.btn, 'Please input data!')
    )
    spe.plt()  
    
  })
  
  # Display node and edge counts below network plot -------------------------
  output$netstats <- renderText({
    validate(
      need(input$load.btn, 'Please input data!')
    )
    paste("Nodes:", gorder(net.obj()), "<br>", "Edges:", gsize(net.obj()))
    
  })
  
  # Recreate taxonomy for barplot -------------------------------------------
  tax.obj.2 <- reactive({
    get.tax(
      tax = tax_data(),
      rank = input$tax.summary,
      graph = net.obj(),
      matrix = use.otu(),
      plot.selection = "Taxonomy"
    )
  })
  
  # Recreate env correlation for boxplot ------------------------------------
  env.cor.2 <- reactive({
    get.env.cor(
      meta = use.env(),
      use.env = input$env.summary,
      matrix = use.otu(),
      plot.selection = "Environmental parameter"
    )
  })
  
  # Add taxonomy and correlation plot information to tidy graph -------------

  net.tbl.2 <- reactive({
    if(!is.null(env.cor.2())){
      tidy.net(
        net = net.obj(),
        matrix = use.otu(),
        clust = clust.obj(),
        env = env.cor.2(),
        tax = tax.obj.2()
      )
    }
    else{
      tidy.net(
        net = net.obj(),
        matrix = use.otu(),
        clust = clust.obj(),
        env = NULL,
        tax = tax.obj.2()
      )
    }
  })
  
  # Taxonomy barplot --------------------------------------------------------
  tax.bar <- reactive({
    tax.plot(
      graph = net.tbl.2(),
      option = input$sum.theme,
      yAxisTitle = input$tax.axis,
      legendTitle = input$tax.summary,
      yAxisData = input$tax.axis,
      env = isolate(use.env()),
      env.var = isolate(input$env.network),
      env.type = env.col(),
      plot.selection = input$colors,
      env.main = env.cor(),
      use.relab = use.otu(),
      use.tax = use.tax()
    )
  })
  
  # Environmental correlation boxplot ---------------------------------------
  cor.obj <- reactive({
    cor.plot(
      graph = net.tbl.2(),
      parameter = isolate(input$env.summary),
      option = input$sum.theme,
      Correlation = isolate(env.cor.2()),
      clus = isolate(input$nclust),
      env = isolate(use.env()),
      env.var = isolate(input$env.network),
      env.type = env.col(),
      plot.selection = input$colors,
      env.main = env.cor()
    )
  })
  
  # Taxonomy and correlation subplot ----------------------------------------
  output$sumplot <- suppressWarnings(renderPlotly({
    if(!is.null(input$env.file)|input$indata == "MGnify"| (input$indata == "Phyloseq" & !is.null(env.cor.2()))){
      if(!(env.col()=="max"&input$colors=="Environmental parameter")&!input$env.summary=="No numeric variables"){
        if(input$colors == "Cluster - kmeans"){
        subplot(style(tax.bar(), showlegend = TRUE),
                style(cor.obj(), showlegend = FALSE),
                margin = 0.1,
                titleX = TRUE,
                titleY = TRUE,
                nrows = 1) %>%
          layout(showlegend=T,height = 800,width = 650+(50*input$nclust) ,legend = list(x = 0.5, y = -0.28,xanchor = "center",orientation = "h",title = "")) %>%
           config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d","zoom", "pan", "select", "zoomIn", "zoomOut","autoScale", "resetScale","lasso2d"))
          } else{
          subplot(style(tax.bar(), showlegend = TRUE),
                  style(cor.obj(), showlegend = FALSE),
                  margin = 0.1,
                  titleX = TRUE,
                  titleY = TRUE,
                  nrows = 1) %>%
            layout(showlegend=T,height = 800,width = 650,legend = list(x = 0.5, y = -0.28,xanchor = "center",orientation = "h",title = "")) %>%
            config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d","zoom", "pan", "select", "zoomIn", "zoomOut","autoScale", "resetScale","lasso2d"))
          }
      } else{
        tax.bar() %>%
        layout(showlegend=T,height = 800+200*min(max(nchar(as.character(env_data()[[input$env.network]])),na.rm = T)/15,2),width = 650+50*length(unique(env_data()[[input$env.network]])),legend = list(x = 0.5, y = -0.1+ max(-max(nchar(as.character(env_data()[[input$env.network]])),na.rm = T)/40,-2),xanchor = "center",orientation = "h",title = "")) %>%
        config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d","zoom", "pan", "select", "zoomIn", "zoomOut","autoScale", "resetScale","lasso2d"))
        }
      } else{
        validate(
          need(input$load.btn, 'Please input data!')
        )
        if(env.col()=="max"){
          tax.bar() %>%
            layout(showlegend=T,height = 800+200*min(max(nchar(as.character(env_data()[[input$env.network]])),na.rm = T)/15,2),width = 650+50*length(unique(env_data()[[input$env.network]])),legend = list(x = 0.5, y = -0.1 + max(-max(nchar(as.character(env_data()[[input$env.network]])),na.rm = T)/40,-2),xanchor = "center",orientation = "h",title = "")) %>%
            config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d","zoom", "pan", "select", "zoomIn", "zoomOut","autoScale", "resetScale","lasso2d"))
        } else{
          tax.bar() %>%
            layout(showlegend=T,height = 800,width = 300,legend = list(x = 0.5, y = -0.2,xanchor = "center",orientation = "h",title = "")) %>%
            config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d","zoom", "pan", "select", "zoomIn", "zoomOut","autoScale", "resetScale","lasso2d"))
        }
   } 
  }))
  

  # Download raw files ------------------------------------------------------
  # Correlation matrix download
  output$dl.matrix <- downloadHandler(
    filename = function() {
      paste("Correlation_matrix_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.csv(cor.mat(), file, quote = FALSE)
    },
    contentType = "text/csv"
  )
  
  # Node data download

  output$dl.aesthetic <- downloadHandler(
    filename = function() {
      paste("Node_data_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      net.tbl()%>%
        activate(nodes) %>%
        as_tibble() %>%
        write.table(file, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
    },
    contentType = "text/csv"
  )
  
  # Edge data download
  output$dl.network <- downloadHandler(
    filename = function() {
      paste("Edge_data_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      net.tbl()%>%
        activate(edges) %>%
        as_tibble() %>% 
        write.table(file, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
    },
    contentType = "text/csv"
  )
  
  # View preview of raw data files ------------------------------------------
  output$rawtables <- renderDataTable({
    validate(
      need(input$load.btn, 'Please input data!')
    )
    data_preview <- reactive({
      if (input$preview == "Node data") {
        node_data <- net.tbl() %>%
          activate(nodes) %>%
          as_tibble()
        return(node_data[1:15,])
      }
      if (input$preview == "Edge data") {
        edge_data <- net.tbl() %>%
          activate(edges) %>%
          mutate(to_name = .N()$name[to], 
                 from_name = .N()$name[from]) %>% 
          as_tibble() %>% 
          select(from = from_name, to = to_name, weight)
        return(edge_data[1:15,])
      }
      if (input$preview == "Similarity matrix") {
        return(cor.mat()[1:15,])
      }
    })
    data_preview()
  }, options = list(scrollX = TRUE))
  
}

shinyApp(ui = ui, server = server)


