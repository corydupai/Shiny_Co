interactive_graph <- function(full_net, subnet_name, label_view = "centroid"){
  
  if(label_view == "centroid"){
    subnet <- full_net %>%
      activate(nodes) %>%
      mutate(id = name, 
             title = gene,
             label = name)
  } else {
    subnet <- full_net %>%
      activate(nodes) %>%
      mutate(id = name, 
             title = name,
             label = gene)
  }
  
  if(grepl("centroid",subnet_name)){
    print("testing")
    subnet_name <- subnet %>%
      activate(nodes) %>%
      as_tibble() %>%
      filter(name == subnet_name) %>%
      select(network) %>%
      unlist() %>%
      unname() %>%
      unique()
    print(subnet_name)
  }
  
  subnet <- subnet %>%
    activate(nodes) %>%
    filter(network == subnet_name) %>%
    mutate(node_filter = centrality_degree()) %>%
    filter(node_filter > 1) %>%
  activate(edges) %>%
  mutate(from = from_node,
         to = to_node)
  
  
vn <- subnet %>% activate(nodes) %>% as_tibble() %>%
  mutate(
    physics = TRUE)
ve <- subnet %>% activate(edges) %>% as_tibble()%>%
  mutate(
    from = from_node,
    to = to_node,
    # color = "black",
    # width = 2,
    physics = FALSE)
visNetwork(vn,edges = ve) %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, collapse = TRUE) %>%
  visInteraction(navigationButtons = TRUE)
}


enrich_plots_func <- function(enriched_dt, term_type){
  print(term_type)
  enrich_plot <- enriched_dt %>% 
    filter(Domain == term_type & qvalue != 0) %>% 
    ggplot(aes(y = ID, x = network, fill = `-log(qvalue)`)) +
    geom_tile(color = "black", size = 1) +
    # scale_fill_viridis(limits = c(1,25))+
    scale_fill_viridis()+
    scale_x_discrete(expand = c(0,0),
                     position = "top",
                     name = "Sub Network\n") +
    scale_y_discrete(expand = c(0,0))+
    theme_cowplot(12) +
    theme(panel.border = element_rect(colour = "black",
                                      size = 0.5),
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          axis.title = element_blank(),
          title = element_text(hjust = 0)
    )
  # plotPNG(width = 400, height = 500, type = "png", file = "/stor/work/Wilke/tingle/WGCNA_expanded/K.pneumoniae/red.png")
}