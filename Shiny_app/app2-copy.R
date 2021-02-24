library(shiny)
library(shinyFiles)
library(shinythemes)
# library(shinydashboard)
library(tidyverse)
library(ggplot2)
library(ggraph)
library(tidygraph)

library(fs)
library(rhandsontable)
library(shinyjs)
library(cowplot)
library(ggrepel)
library(DT)
# library(cairoDevice)

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
                    choices = list("choose",
                                   "Bacillus subtilis",
                                   "Escherichia coli",
                                   "Klebsiella pneumonia",
                                   "Pseudomonas aeruginosa"
                                   )),  # End of selectizeInput()
      
            hidden(
              div(id = "choose_subnetwork",
                  selectizeInput("subnetwork",
                                 label = "Subnetwork Selection",
                                 choices=NULL))), # End of hidden()
          
            hidden(
              div(id = "choose_graph_o_enrich",
                  selectizeInput("graph_o_enrich",
                                 label = "Select what you would like to visualize",
                                 choices=NULL))),
              hidden(
                div(id = "submit_Button",
                    column(12,
                           h4("Submit")))) # End hidden 4

      ), # end of sidebarPanel()
    
####### Main panel #######
      mainPanel(
        
        
        # hidden(
        #   div(id = "graph_plot",
        #       column(12,
        #              h4("Here's your graph"),
        #              plotOutput("graph")))), # End hidden 1
        # 
        # hidden(
        #   div(id = "graph_table",
        #       column(12,
        #              h4("Network loci"),
        #              dataTableOutput("gtable")))), # End hidden 2
        
        hidden(
          div(id = "graph_tabs",
              column(12,
                     h4(""),
                     tabsetPanel(type = "tabs",
                                 tabPanel("Graph Subnetwork", plotOutput("graph")),
                                 tabPanel("Network Loci", dataTableOutput("gtable"))
                                 # tabPanel("Table", tableOutput("table"))
                     )))),
        
        
        # hidden(
        #   div(id = "enrichment_plot",
        #       column(12,
        #              h4("Here's your plot"),
        #              imageOutput("enrichment")))), # End hidden 3
        # 
        # 
        # hidden(
        #   div(id = "annot_table",
        #       column(12,
        #              h4("KEGG and GO annotations"),
        #              dataTableOutput("annotations")))), # End hidden 4
        
        hidden(
          div(id = "enrich_tabs",
              column(12,
                     h4(""),
        tabsetPanel(type = "tabs",
                    tabPanel("Enrichment Plot", imageOutput("enrichment")),
                    tabPanel("Annotations Table", dataTableOutput("annotations"))
                    # tabPanel("Table", tableOutput("table"))
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
  # 
  updateSelectizeInput(session, 'subnetwork', choices = c("choose", "all"), server = TRUE)
  updateSelectizeInput(session, 'graph_o_enrich', choices = c("choose","Graph", "Gene Enrichment"), server = TRUE)
  
  # path_tr <- eventReactive(input$maindir, {
  #   if(is.numeric(input$maindir)) return(NULL)
  #   val_out <- "/stor/work/Wilke/tingle/WGCNA_expanded/K.pneumoniae"
  #   val_out
  # })
  
species_tr <- eventReactive(input$choose_species,{
    if(input$choose_species == "choose") return(NULL)
    if(input$choose_species == "Klebsiella pneumonia")
    {
      val_out <- "/stor/work/Wilke/tingle/WGCNA_expanded/K.pneumoniae/"
      # show("choose_subnetwork")
      return(val_out)
      
    }
    return(NULL)
       # return(TRUE)
  })
  
subnetwork_tr <- eventReactive(input$subnetwork,{
  if(input$subnetwork == "all")
    return(dynamicColors)
  else
    return(as.character(input$subnetwork))
  
})

graph_o_enrich_tr <- eventReactive(input$graph_o_enrich,{
    return(input$graph_o_enrich)
})
  
  observeEvent(species_tr(),{
    print(species_tr())
    show("choose_subnetwork")
    
    observeEvent(subnetwork_tr(),{
      if(subnetwork_tr() == as.character(input$subnetwork))
      {
        show("choose_graph_o_enrich")
      }
    })
    
    observeEvent(graph_o_enrich_tr(),{
      if(graph_o_enrich_tr() == "Graph")
      {
        show("graph_tabs")
        # show("graph_plot")
        # show("graph_table")
      }
    })
    
    observeEvent(graph_o_enrich_tr(),{
      if(graph_o_enrich_tr() == "Gene Enrichment")
      {
        show("enrich_tabs")
        # show("enrichment_plot")
        # show("annot_table")
        # show("graph_table")
      }
    })
    
    
    }
  )

  full_graph <- eventReactive(species_tr(),{
    f_graph <- readRDS(paste0(species_tr(), "/filtered_graph_0.10.RDS")) %>%
      activate(nodes) %>%
      filter(dynamicColors == subnetwork_tr())
   
  } #end of eventReactive {f_graph}
  ) #end of full_graph
  
  output$graph <- 
    renderPlot({
      req(full_graph())
      
      if(length(subnetwork_tr()) == 1){
        ggraph(full_graph(), layout="stress") +
          geom_edge_link(alpha = 0.1) + 
          geom_node_point(aes(color=dynamicColors), 
                          colour="black",
                          fill = subnetwork_tr(), 
                          show.legend=FALSE,
                          size = 4,
                          shape = 21) +
          geom_node_text(aes(label = name), size=4, repel = TRUE)
      }
      
    }) # End of graph
  
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

