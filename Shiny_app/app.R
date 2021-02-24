library(shiny)
library(shinyFiles)
library(shinythemes)
library(tidyverse)
library(ggplot2)
library(tidygraph)
library(fs)
library(rhandsontable)
library(shinyjs)
library(cowplot)
library(ggrepel)
library(DT)
# install.packages("visNetwork")
library(visNetwork)
library(data.table)
library(plotly)
library(viridis)
source("Rcode/functions.R")
files_dir <- "/stor/work/Wilke/tingle/WGCNA_expanded/bacteria/"
theme_set(theme_cowplot(12))
Sys.umask("000")

########################################################
# Define the UI
########################################################
ui <- 
  shinyUI(fluidPage(
    
    theme = shinytheme("sandstone"),
    
    titlePanel("WGCN Visualization"),
    
    sidebarLayout(position = "left", 
                  
                  ####### Sidebar #######
                  sidebarPanel(
                    
                    helpText("Visualize Weighted Gene Co-Expression Networks and Annotations"),
                    
                    useShinyjs(), 
                    
                    selectizeInput("choose_species", 
                                   label = "Select a species of bacteria",
                                   choices = NULL),  # End of selectizeInput()
                    hidden(
                      div( id = "choose_chooser",
                           radioButtons("chooser",
                                        label = "How would you like to select your subnetwork?",
                                        choices = c("Name",
                                                    "Centroid"),
                                        selected = character(0),
                                        inline = TRUE)
                      )
                    ),
                    hidden(
                      div(
                        id = "choose_subnetwork",
                        selectizeInput("subnetwork",
                                       label = NULL,
                                       choices = NULL
                        )
                      )
                    ),
                    
                    hidden(
                      div(id = "submit_Button",
                          column(12,
                                 h4("Submit")))) # End hidden 4
                    
                  ), # end of sidebarPanel()
                  
                  ####### Main panel #######
                  mainPanel(
                    id = "MAIN",
                    column(12,
                           h4(""),
                           hidden(
                             div( 
                               id = "graph_tabs",
                               tabsetPanel(
                                 type = "tabs",
                                 tabPanel("Gene Mapping",
                                          hidden(
                                            div( id = "Strains",
                                                 selectizeInput("strain_select",
                                                                label = "Select your strain of interest",
                                                                choices = NULL))),
                                          dataTableOutput("gene_info")
                                 ),
                                 tabPanel("Enriched Terms",
                                          tabsetPanel(type = "tabs",
                                                      tabPanel("KEGG Pathway",
                                                               plotlyOutput("KEGG_terms", height = "800px")),
                                                      tabPanel("KEGG Enzyme",
                                                               plotlyOutput("Enzyme_terms", height = "800px")),
                                                      tabPanel("Biological Process",
                                                               plotlyOutput("BP_terms", height = "800px")),
                                                      tabPanel("Molecular Function",
                                                               plotlyOutput("MF_terms", height = "800px")),
                                                      tabPanel("Cellular Component",
                                                               plotlyOutput("CC_terms", height = "800px"))
                                                      )
                                          
                                 ),
                                 
                                 
                                 tabPanel("Network Visualization", 
                                          visNetworkOutput("graph"),
                                          selectizeInput("select_labels",
                                                         label = "Select how you want the graph labelled",
                                                         choices=c("centroid",
                                                                   "gene ID"),
                                                         selected = "centroid")
                                          
                                 ),
                                 tabPanel("Coexpression By Gene", 
                                          selectizeInput("COI",
                                                         label = "Centroid of interest",
                                                         choices = NULL),
                                          dataTableOutput("coexpressed_genes")),
                                 id = "Coexp_tabs"
                               ))
                           )),
                    
                    hidden(
                      div(id = "enrich_tabs",
                          column(12,
                                 h4(""),
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Enrichment Plot", imageOutput("enrichment")),
                                             tabPanel("Annotations Table", dataTableOutput("annotations"))
                                 )))) # End hidden 5
                    
                    
                    
                  ) # end of mainPanel()
                  
    ) # send sidebarLayout()
  ) # end fluidPage()
  ) # end shinyUI()



########################################################
# Define the Server
########################################################

server <- function(input, output, session){
  # volumes <- c(Home = "/stor/work/Wilke/tingle/WGCNA_expanded")
  # shinyFileChoose(input, "DE", roots = volumes, session = session)
  # shinyFileChoose(input, "fastqs_exist", roots = volumes, session = session)
  # shinyFileChoose(input, "tr_exists", roots = volumes, session = session)
  # shinyDirChoose(input, "maindir", roots = volumes, session = session, restrictions = system.file(package = "base"))
  hideTab(inputId = "Coexp_tabs", target = "Network Visualization")
  hideTab(inputId = "Coexp_tabs", target = "Gene Mapping")
  hideTab(inputId = "Coexp_tabs", target = "Coexpression By Gene")
  species_list <- c(
    " " = "choose",
    "Bacillus subtilis" = "B.subtilis",
    "Escherichia coli" = "E.coli",
    "Klebsiella pneumonia" = "K.pneumoniae",
    "Pseudomonas aeruginosa"  = "P.aeruginosa"
  )
  updateSelectizeInput(session, 'choose_species', choices = species_list, server = TRUE)
  updateSelectizeInput(session, 'graph_o_enrich', choices = c("choose","Graph", "Gene Enrichment"), server = TRUE)
  big_graph <- eventReactive(input$choose_species,{
    if(input$choose_species == "" | input$choose_species == "choose") {
      return(NULL)
    } else {
      strains <- read_tsv(paste0(files_dir,input$choose_species,"/data4app/genome_list.csv"),
                          col_names = FALSE) %>%
        unlist() %>%
        unname()
      
      updateSelectizeInput(session, 'strain_select', choices = strains, server = TRUE)
      show("Strains")
      graph_tbl <- readRDS(paste0(files_dir,input$choose_species,"/data4app/filtered_graph_0.10.RDS"))
      return(graph_tbl)
    }
  })
  
  
  observeEvent(big_graph(),{
    # print(big_graph())
    subnets <- big_graph() %>%
      activate(nodes) %>%
      as_tibble() %>%
      select(network,dynamicColors) %>%
      unique()
    
    # updateSelectizeInput(session, 'subnetwork', choices = subnets, server = TRUE)
    show("graph_tabs")
    showTab(inputId = "Coexp_tabs", target = "Gene Mapping")
    show("choose_chooser")
    enrich_dt <- fread(paste0(files_dir,input$choose_species,"/data4app/enriched_terms.csv")) %>%
      left_join(subnets)
    
    
    output$KEGG_terms <- renderPlotly(
      enrich_plots_func(enrich_dt , "KEGG Pathways") %>% ggplotly() )
    output$Enzyme_terms <- renderPlotly(
      enrich_plots_func(enrich_dt , "KEGG Enzymes") %>% ggplotly() )
    output$CC_terms <- renderPlotly(
      enrich_plots_func(enrich_dt , "Cellular Component") %>% ggplotly() )
    output$MF_terms <- renderPlotly(
      enrich_plots_func(enrich_dt , "Molecular Function") %>% ggplotly() )
    output$BP_terms <- renderPlotly(
      enrich_plots_func(enrich_dt , "Biological Process") %>% ggplotly() )
  })
  
  
  
  observeEvent(c(input$chooser),{
    req(big_graph())
    COIs <- big_graph() %>%
      activate(nodes) %>%
      as_tibble()
    
    SNs <- COIs %>%
      select(network) %>%
      arrange(network) %>%
      unique() %>%
      unlist() %>%
      unname()
    
    COIs <- COIs %>%
      select(name) %>%
      unlist() %>%
      unname()
    
    updateSelectizeInput(session, 'COI', choices = COIs, server = TRUE)
    show("choose_subnetwork")
    if(input$chooser == "Name"){
      updateSelectizeInput(session, 'subnetwork', choices = SNs, server = TRUE, label = "Network Name")
    } else {
      updateSelectizeInput(session, 'subnetwork', choices = COIs, server = TRUE, label = "Centroid ID")
    }
    
    
  })
  
  graph_o_enrich_tr <- eventReactive(input$graph_o_enrich,{
    return(input$graph_o_enrich)
  })
  
  observeEvent(
    c(input$subnetwork,
      input$select_labels),{
        req(input$choose_species)
        show("choose_gene")
        show("choose_graph_o_enrich")
        if(grepl("centroid", input$subnetwork)){
          
          GOIs <- big_graph() %>%
            activate(nodes) %>%
            filter(name == input$subnetwork) %>%
            as_tibble() %>%
            select(network) %>%
            unlist() %>%
            unname()
          
          # print(GOIs)
          COIs <- big_graph() %>%
            activate(nodes) %>%
            as_tibble() %>%
            filter(network == GOIs) %>%
            select(name) %>%
            unlist() %>%
            unname() %>%
            unique()
          # print(COIs)
          updateSelectizeInput(session, 'COI', choices = COIs, selected = input$subnetwork)
        }
        output$graph <-
          renderVisNetwork(interactive_graph(big_graph(),input$subnetwork, input$select_labels))
        showTab(inputId = "Coexp_tabs", target = "Network Visualization")
        showTab(inputId = "Coexp_tabs", target = "Coexpression By Gene")
        
      })
  
  
  observeEvent(
    input$COI,{
      if(input$COI != ""){
        temp_vals <-
          big_graph() %>%
          activate(edges) %>%
          filter(from_node == input$COI |
                   to_node == input$COI) %>%
          as_tibble() %>%
          mutate(to_node =
                   if_else(to_node == input$COI,
                           from_node,
                           to_node),
                 from_node = input$COI,
                 coexpression = weight) %>%
          select(from_node,
                 to_node,
                 coexpression) %>%
          left_join(big_graph() %>%
                      activate(nodes) %>%
                      as_tibble() %>%
                      rename("to_gene" = gene) %>%
                      select(name,to_gene,network),
                    by = c("to_node" = "name")) %>%
          left_join(big_graph() %>%
                      activate(nodes) %>%
                      as_tibble() %>%
                      rename("from_gene" = gene) %>%
                      select(name,from_gene),
                    by = c("from_node" = "name"))
        output$coexpressed_genes <- renderDataTable(
          
          DT::datatable(temp_vals,
                        filter = "top",
                        extensions = 'Buttons',
                        options = list(
                          dom = 'Bfrtip',
                          
                          buttons = c('csv', I('colvis')),
                          scrollX = TRUE,
                          fixedColumns = list(leftColumns = 2, rightColumns = 1),
                          pageLength = 10,
                          lengthMenu = c(5, 10, 25, 50)
                        )
          )
        )
      }
    }
  )
  
  output$gtable <- renderDataTable(
    DT::datatable(full_graph() %>% 
                    activate(nodes) %>%
                    as_tibble(),
                  filter = "top",
                  extensions = 'Buttons',
                  options = list(
                    dom = 'Bfrtip',
                    
                    buttons = c('csv', I('colvis')),
                    scrollX = TRUE,
                    fixedColumns = list(leftColumns = 2, rightColumns = 1),
                    pageLength = 5,
                    lengthMenu = c(5, 10, 25, 50)
                  )
    )
  ) # End of gtable
  
  output$gene_info <- renderDataTable(
    DT::datatable(
      fread(paste0(files_dir,input$choose_species,"/data4app/locus_table.csv"),
                        select = c("centroid", input$strain_select)) %>%
        right_join(
          big_graph() %>%
                    activate(nodes) %>%
                    as_tibble() %>%
      select(name, "gene name" = gene, network),
        by = c("centroid" = "name")
        ),
    filter = "top",
    extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      
      buttons = c('csv', I('colvis')),
      scrollX = TRUE,
      pageLength = 10,
      lengthMenu = c(5, 10, 25, 50)
    )
  )
  )
  
  output$enrichment <- 
    renderImage({
      
      # tag$img(src = "/stor/work/Wilke/tingle/WGCNA_expanded/K.pneumoniae/red.png",
      #          alt = "Placeholder enrichment plot",
      #          height = "75%", 
      #          align="left")
      list(src = "/stor/work/Wilke/tingle/WGCNA_expanded/K.pneumoniae/darkred.png",
           # contentType = 'png',
           width = 600,
           height = 800,
           # height = "100%",
           alt = "This is alternate text")
      
    }, deleteFile = FALSE) # End of graph
  
  
  output$annotations <- renderDataTable(
    DT::datatable(read_csv(paste0(species_tr(), "/test_annot.csv")) %>% 
                    filter(dynamicColors == subnetwork_tr()) %>% 
                    select(Description, geneID, Domain, qvalue, 
                           `-log(qvalue)`, Count, full_counts, BR), 
                  # activate(nodes) %>%
                  # as_tibble(),
                  filter = "top",
                  extensions = 'Buttons',
                  options = list(
                    dom = 'Bfrtip',
                    
                    buttons = c('csv', I('colvis')),
                    scrollX = TRUE,
                    fixedColumns = list(leftColumns = 2, rightColumns = 1),
                    pageLength = 5,
                    lengthMenu = c(5, 10, 25, 50)
                  )
    )
  ) # End of annotations
  
  
 
  
  
  
  
} # End of server function

# Run the application 
shinyApp(ui = ui, server = server)

