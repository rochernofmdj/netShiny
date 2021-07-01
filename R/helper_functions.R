library(scales)
#Function that gives you the ability to manipulate the width and color of
#the edges in an igraph network
#Credit: Gábor Csárdi
#https://lists.gnu.org/archive/html/igraph-help/2013-03/msg00030.html
mycircle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x = x, y = y, bg = bg, fg = fg, lwd = lwd,
                   circles = size, add = TRUE, inches = FALSE)
         })
}

igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                         plot = mycircle, parameters = list(vertex.frame.color = 1, vertex.frame.width = 1))

# Not in operator
`%nin%` <- Negate(`%in%`)

# Arguments of netphenogeno and selectnet
args_netphenogeno <- formalArgs(netgwas::netphenogeno)[-1]
args_netphenogeno <- args_netphenogeno[-length(args_netphenogeno)]
args_selectnet <- formalArgs(netgwas::selectnet)[-1]
args_selectnet <- args_selectnet[-length(args_selectnet)]

# Custom colors that we (semi) manually selected to use in the plots
# These colors were chosen such that any two colors next to each are as distinct as possible
custom.col <- c("#DBB165",  "#52854C", "#4E84C4", "#C3D7A4","#C4961A", "#8B4513", "#293352",
                "#E69F00", "#56B4E9", "#009E73", "#B62A3D", "#0072B2", "#D55E00", "#CC79A7",
                "#000000", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
                "#B15928", "#FFD700", "#C5B358")

#Custom bsModal that does not close when users click outside the dialog modal
#or when users press the escape key
custombsModal <- function(...) {
  bmodal <- shinyBS::bsModal(...)
  bmodal[[2]]$`data-backdrop` <- "static"
  bmodal[[2]]$`data-keyboard` <- "false"
  return(bmodal)
}

#This function checks if only a few nodes are chosen to create a subnetwork
#If a subnetwork needs to be created, it does so by only getting the nodes that are chosen
#with their neighbors
get_subnet <- function(mat){
  if(!is.null(vals$subgraph_nodes)){
    ind_subgraph_nodes <- which(dimnames(mat)[[2]] %in% vals$subgraph_nodes)
    summ_df <- summary(mat)
    summ_df <- summ_df[summ_df$i %in% ind_subgraph_nodes | summ_df$j %in% ind_subgraph_nodes, ]
    nms_subgraph_nodes <- sort(union(summ_df$i, summ_df$j))
    nms_subgraph_nodes <- dimnames(mat)[[2]][nms_subgraph_nodes]
    mat <- mat[nms_subgraph_nodes, nms_subgraph_nodes]
  }
  return(mat)
}

#Function that returns an adjacency or weighted matrix with only nodes
#that satisfy threshold value(s) given. If (both) threshold value(s) are NULL, it will
#return an adjacency or weight matrix that only include non-isolated nodes
getNZ <- function(mat, ntraits = NULL, th_t = NULL, th_m = NULL){
  shiny::validate(shiny::need(!is.null(mat), "Getting proper data"))
  mat <- Matrix::Matrix(mat)
  n_nodes <- length(vals$node_names)

  if(is.null(ntraits)){

    if(!is.null(th_t)){
      mat[abs(mat) < th_t] <- 0
    }
    diag(mat) <- 0 #Makes diagonal values 0
    mat <- Matrix::drop0(mat, tol = 0, is.Csparse = TRUE) #Returns sparse matrix with no explicit zeroes (including removing diagonal zeroes)
    nzvec <- which(colSums(mat != 0) == 0) #Gets indices of isolated nodes
    mat <- mat[-nzvec, -nzvec]
    invisible(mat)
  }

  else{

    if (!is.null(th_m)){
      mat[(ntraits+1):n_nodes, (ntraits+1):n_nodes][abs(mat[(ntraits+1):n_nodes, (ntraits+1):n_nodes]) < th_m] <- 0
    }

    if (!is.null(th_t)){
      mat[, 1:ntraits][abs(mat[, 1:ntraits]) < th_t] <- 0
      mat[1:ntraits, ][abs(mat[1:ntraits, ]) < th_t] <- 0
    }

    diag(mat) <- 0  #Makes diagonal values 0

    mat <- Matrix::drop0(mat, tol = 0, is.Csparse = TRUE) #Returns sparse matrix with no explicit zeroes for opt.adj, thus removing diagonal 0 values
    nzvec <- which(colSums(mat != 0) == 0) #gets indices of traits and markers that do not have any links
    nzvec <- nzvec[nzvec > ntraits] #gets indices of only markers that do not have any links
    mat <- mat[-nzvec, -nzvec]
    invisible(mat)
  }
}

#This function takes as argument an adjacency or weight matrix
#The function returns the names of only non-isolated nodes
get_nz_nodes <- function(mat, ntraits = NULL){
  mat <- getNZ(mat, ntraits)
  nz_names <- dimnames(mat)[[2]]
  return(nz_names)
}

#Function that get the needed igraph object to use to make a network
#This function takes in a (sparse) matrix and returns an igraph object
get_g_complete <- function(mat){
  g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
  #In gxe mode traits and marker nodes have different sizes
  if(vals$mode == "gxe" && !is.null(vals$map_nodes)){
    labs <- vals$map_nodes[vals$map_nodes$node %in% igraph::vertex_attr(g, "name"), ]
    traits_mat <- sum(labs$node_group == "Trait")
    igraph::V(g)$size <- c(rep(20, traits_mat), rep(14, length(labs$node_group) - traits_mat))
  }

  #Assign color to edges depending if it is negative (red) or (positive) blue
  if(!is.null(igraph::E(g)$weight)){
    edge_cols <- ifelse(igraph::E(g)$weight > 0, "blue", "red")
    value <- abs(igraph::E(g)$weight)
    igraph::E(g)$color <- edge_cols
    igraph::E(g)$width <- value
  }

  return(g)
}

#Function that gets the proper layour for the networks
#This function takes a (sparse) matrix as an argument and
#returns a dataframe with the coordinates layout
get_igraph_lay <- function(mat){
  #Get an initial set of coordinates that later will be matched
  #with the coordinates of the global network
  mat@x[abs(mat@x) > 0] <- 1
  g_lay <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected")
  lay <- igraph::layout_nicely(g_lay)
  row.names(lay) <- names(igraph::V(g_lay))

  #Changing coordinates to match global network
  for (i in 1:nrow(lay)){
    if(is.element(rownames(lay)[i], rownames(coords()))){
      lay[rownames(lay)[i], ] <- coords()[rownames(lay)[i], ]
    }
  }

  #Make sure that isolated nodes are not too far away from the rest of the nodes
  if (shiny::isolate(input$layout == "Tree")){
    s <- which(igraph::degree(g_lay) == 0)
    s <- names(s)
    lay_y <- sort(coords()[, 2])
    new_y <- lay_y[length(s) + 1]
    new_xs <- coords()[, 1][coords()[, 2] == new_y]
    new_X <- max(new_xs)
    for (i in s){
      new_X <- new_X + 3
      lay[i, ][1] <- new_X
      lay[i, ][2] <- new_y
    }
  }

  return(lay)
}

#This function gets a visnetwork to be shown in the app
#This function takes in a (sparse) matrix, igraph object, and layout dataframe
#and returns a visnetwork
get_vis_net <- function(mat, g, lay){
  sel_nodes <- dimnames(mat)[[1]]
  test.visn <- visNetwork::toVisNetworkData(g)
  sel_by <- NULL

  if(!is.null(vals$map_nodes)){
    if(vals$mode == "gxe"){
      sel_by <- "Chromosome"
      test.visn$nodes$Chromosome <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_group
    }
    else{
      sel_by <- "Group"
      test.visn$nodes$Group <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_group
    }
    test.visn$nodes$color.background <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_color
    igraph::V(g)$color <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_color
  }

  if(!is.null(igraph::E(g)$weight)){
    test.visn$edges$value <- abs(igraph::E(g)$weight)
  }
  test.visn$nodes$font.size <- 17

  #vals$g1 <- g

  visNetwork::visNetwork(test.visn$nodes, test.visn$edges, height = '1000px', width = '1000px') %>%
    visNetwork::visIgraphLayout(layout = "layout.norm", layoutMatrix = lay, randomSeed = vals$rseed) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = T, hover = T), selectedBy = sel_by) %>%
    visNetwork::visInteraction(multiselect = TRUE) %>%
    visNetwork::visEvents(doubleClick =  "function(nodes){Shiny.onInputChange('click', nodes.nodes[0]);;}") %>%
    visNetwork::visEdges(smooth = list(enabled = TRUE, type = "curvedCCW", roundness = input$roundness)) %>%
    return
}

get_con <- function(mat, sett, node){
  dimnms <- dimnames(mat)[[2]]
  node_ind <- which(dimnms == node) #Get index of chosen node
  summ <- summary(mat)
  summ <- summ[summ$i == node_ind | summ$j == node_ind, ]
  summ <- summ[!(summ$i == node_ind & summ$j == node_ind), ]
  summ$setting <- rep(sett, nrow(summ))
  summ$cons <- ifelse(summ$i == node_ind, dimnms[summ$j], dimnms[summ$i])
  summ <- summ[c("cons", "x", "setting")]
  colnames(summ) <- c("cons", "value", "setting")
  return(summ)
}

#This function assigns color and label to each node
assign_col_lab <- function(type_obj, lab_or_col, n_nodes){
  if (lab_or_col == "col"){
    ln <- length(type_obj)
    cols <- mapply(rep, custom.col[1:ln], unlist(type_obj), SIMPLIFY = FALSE)
    cols <- unlist(cols)
    names(cols) <- NULL
    if(!is.null(vals$ntraits)){
      cols <- c(cols, rep(custom.col[ln+1], n_nodes))
    }

    return(cols)
  }
  else{
    labs <- mapply(rep, names(type_obj), unlist(type_obj), SIMPLIFY = FALSE)
    labs <- unlist(labs)
    names(labs) <- NULL
    if(!isnull(vals$ntraits)){
      labs <- c(labs, rep("Markers", n_nodes))
    }
    return(labs)
  }
}

mapping_to_string <- function(df){
  #If less unique values in first column, then first column contains groups for nodes, or vice versa
  bool_group_col <- length(unique(df[, 1])) < length(unique(df[, 2]))
  #A list of the groups with the nodes corresponding to each group
  if(bool_group_col){
    groups <-  split(x = df[, 2], f = df[, 1])
  }
  else{
    groups <-  split(x = df[, 1], f = df[, 2])
  }
  n_members <- unlist(lapply(groups, length)) #Vector containing number of members for each group
  group_string <- paste(names(groups), n_members, sep = ":", collapse = ", ")
  return(group_string)
}

string_to_df <- function(group_string){
  #Splits the trait types first by comma
  #Then splits them by semicolon, and convert into a dataframe, then back into list
  #So we get a list with the grouping with their frequency as a list
  groups <- trimws(strsplit(group_string, split = ",")[[1]])
  bools <- grepl("\\w+:{1}\\d+", groups) #Checks if strings have proper format
  if (all(bools)){
    groups <- strsplit(groups, ":")
    groups <- do.call(rbind.data.frame, groups)
    colnames(groups) <- c("node_group", "freq")
    groups$freq <- as.numeric(as.character(groups$freq))
    if(sum(groups$freq) != length(all_node_nms)){
      return(0)
    }
    n_groups <- length(groups$node_group)
    node_group <- mapply(rep, groups$node_group, groups$freq, SIMPLIFY = FALSE)
    node_color <- mapply(rep, custom.col[1:n_groups], groups$freq, SIMPLIFY = FALSE)
    df <- data.frame("node" = all_node_nms, "node_group" = unlist(node_group), "node_color" = unlist(node_color))
  }
  return(df)
}

#Function that gets the trait types grouping
#This function takes in a vector and returns a dataframe
#with a trait grouping and its respective total number
get_trait_groups <- function(){
  #Split the trait types first by comma
  #Then split them by semicolon, and convert into a dataframe, then back into list
  #So we get a list with the grouping with their frequency as a list
  trt_typs <- trimws(strsplit(input$trait_types, split = ",")[[1]])
  trt_typs_bool <- gsub("\\s", "", trt_typs)
  bools <- grepl("\\w+:{1}\\d+", trt_typs_bool)
  if (all(bools)){
    trt_typs <- strsplit(trt_typs, ":")
    trt_typs <- do.call(rbind.data.frame, trt_typs)
    colnames(trt_typs) <- c("envs", "freq")
    trt_typs$freq <- as.numeric(as.character(trt_typs$freq))
    trt_typs$envs <- trimws(as.character(trt_typs$envs))
    return(trt_typs)
  }
  else{
    return(NULL)
  }
}

#Function to get the plot for the distribution of the weights of the settings
#This function takes no arguments and returns a ggplot object
get_par_cor_plot <- function(){
  dat <- data.frame(type = factor(),
                    weights = numeric())
  if(vals$mode == "gxe"){
    x_lab <- "Distribution Partial Correlations"
  }
  else{
    x_lab <- "Distribution Edge Weights"
  }
  for (i in 1:length(vals$networks)){
    data <- vals$networks[[i]]
    data <- getNZ(vals$networks[[i]], vals$ntraits, input$cor_t, input$cor_m)
    data <- get_subnet(data)

    shiny::validate(shiny::need(nrow(data) > 0, 'No Connections'))

    diag(data) <- 0
    data <- Matrix::drop0(data, tol = 0, is.Csparse = TRUE)
    lst <- as.vector(data)
    lst <- lst[lst != 0]
    dat <- rbind(dat, data.frame(type = rep(vals$sett_names[i], length(lst)), weights = lst))
  }

  p <- ggplot2::ggplot(data = dat, ggplot2::aes(x = weights, fill = type)) +
    ggplot2::geom_histogram(colour = "black", fill = "#007c00", bins = input$par_cor_bins, ggplot2::aes(y = stat(width*density))) +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab("Proportion") +
    ggplot2::facet_wrap(~type) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = input$par_cor_breaks), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                   axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title.x = ggplot2::element_text(size = 16, face = "bold"),
                   strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "#F5F2F2"),
                   panel.background = ggplot2::element_rect(fill = "#F5F2F2"))
  return(p)
}

#Function that gets the weights analysis plot
#This function takes in no argument and returns a ggplot object
get_weights_analysis_plot <- function(){
  df_conns <- mapply(get_con, vals$networks, vals$sett_names, input$marker, SIMPLIFY = FALSE)
  df_conns <- do.call("rbind", df_conns)
  if(shiny::isTruthy(vals$map_nodes)){
    df_conns <- merge(df_conns, vals$map_nodes[c("node", "node_group")], by.x = "cons", by.y = "node", all.x = TRUE)
    df_conns <- df_conns[order(factor(df_conns$cons, levels = vals$map_nodes$node)),]
    df_conns$cons <- factor(df_conns$cons, levels = vals$map_nodes$node)
    df_conns$node_group <- factor(df_conns$node_group, levels = unique(vals$map_nodes$node_group))
    row.names(vals$map_nodes) <- vals$map_nodes$node
    frame_colors <- unique(vals$map_nodes[unique(as.character(df_conns$cons)), ]$node_color)
  }
  else{
    df_conns$node_group <- df_conns$setting
  }

  lgnd.pos <- ifelse(is.null(vals$map_nodes), "none", "right")
  if(vals$mode == "gxe"){
    node_text <- "Marker: "
    con_text <- "Par. Cor.: "
    group_text <- "Chromosome: "
    y_lab <- "Partial Correlation"
  }
  else{
    node_text <- "Node: "
    con_text <- "Weight: "
    group_text <- "Group: "
    y_lab <- "Weight"
  }

  p <- ggplot2::ggplot(data = df_conns, ggplot2::aes(x = cons, text = paste0(node_text, cons, "\n",
                                                                             con_text, round(value, 3), "\n",
                                                                             group_text, node_group))) +
    ggplot2::geom_point(ggplot2::aes(y = value), color = "red") +
    #geom_line(aes(y = values, group = 1)) +
    ggplot2::geom_bar(ggplot2::aes(weight = value, fill = node_group), show.legend = FALSE) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL) +
    ggplot2::facet_wrap(~setting, nrow = length(vals$networks)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1, size = 10),
                   legend.text = ggplot2::element_text(size = 15)) + #, axis.title.y = element_text(size = 15)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(fill = NULL) +
    ggplot2::theme(legend.position = lgnd.pos)

  if(shiny::isTruthy(vals$map_nodes)){
    p <- p +
      ggplot2::scale_fill_manual(values = frame_colors)
  }

  return(p)
}

getDif <- function(){
  mat1 <- vals$networks[[input$net1]]
  mat2 <- vals$networks[[input$net2]]

  newobj <- list()

  mat1@x[mat1@x < 0] <- -1
  mat1@x[mat1@x > 0] <- 4
  mat2@x[mat2@x < 0] <- -1
  mat2@x[mat2@x > 0] <- 4

  diff <- mat1 - mat2
  diff@x[diff@x %in% c(-5, 5)] <- 1
  diff@x[diff@x %in% c(-1, 4)] <- 2
  diff@x[diff@x %in% c(-4, 1)] <- 3
  diff <- Matrix::drop0(diff)

  diff <- getNZ(diff)
  nms <- dimnames(diff)[[2]]

  clrs <- c("blue", "red", "green")
  summ <- summary(diff)
  summ$i <- nms[summ$i]
  summ$j <- nms[summ$j]
  summ$color <- clrs[summ$x]

  diff@x[diff@x != 0] <- 1
  newobj$opt.adj <- diff
  newobj$changes <- summ

  shiny::validate(
    shiny::need(dim(newobj$opt.adj)[[1]] != 0, "No Difference Between Networks")
  )

  g <- igraph::graph_from_adjacency_matrix(newobj$opt.adj, mode = "undirected")
  igraph::E(g)$color <- newobj$changes$color
  test.visn <- visNetwork::toVisNetworkData(g)

  sel_nodes <- dimnames(newobj$opt.adj)[[1]]
  sel_by <- NULL
  if(!is.null(vals$map_nodes)){
    if(vals$mode == "gxe"){
      sel_by <- "Chromosome"
      test.visn$nodes$Chromosome <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_group
    }
    else{
      sel_by <- "Group"
      test.visn$nodes$Group <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_group
    }
    test.visn$nodes$color.background <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_color
    igraph::V(g)$color <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_color
  }

  ledges <- data.frame(color = c("blue", "green", "red"),
                       label = c("Sign Change", "Gained", "Lost"), arrows = c("unidrected", "undirected", "undirected"))

  visNetwork::visNetwork(test.visn$nodes, test.visn$edges, height = '1000px', width = '1000px') %>%
    visNetwork::visIgraphLayout()  %>%
    visNetwork::visOptions(highlightNearest = list(enabled = T, hover = T), selectedBy = sel_by) %>%
    visNetwork::visLegend(addEdges = ledges, position = "right") %>%
    return
}

get_diff_sets <- function(){
  mat1 <- getNZ(vals$networks[[input$net1]], vals$ntraits, input$cor_t, input$cor_m)
  mat2 <- getNZ(vals$networks[[input$net2]], vals$ntraits, input$cor_t, input$cor_m)

  mat1_nms <- dimnames(mat1)[[2]]
  mat2_nms <- dimnames(mat2)[[2]]

  if(input$sets_selin == "Union"){
    mats_nms <- union(mat1_nms, mat2_nms)
  }

  else if(input$sets_selin == "Intersection"){
    mats_nms <- intersect(mat1_nms, mat2_nms)
  }

  else if(input$sets_selin == "Complement"){
    mats_nms <- setdiff(mat1_nms, mat2_nms)
  }

  shiny::validate(
    shiny::need(length(mats_nms) != 0, "No Difference Between Networks")
  )

  if(shiny::isTruthy(vals$map_nodes)){
    mats_nms_grps <- vals$map_nodes[vals$map_nodes$node %in% mats_nms, ]$node_group
    df_tab <- data.frame(mats_nms, mats_nms_grps)
    if(vals$mode == "gxe"){
      colnames(df_tab) <- c("Markers", "Chromosome")
    }
    else{
      colnames(df_tab) <- c("Node", "Node Group")
    }
  }
  else{
    df_tab <- data.frame("Node" = mats_nms)
  }


  dt <- DT::datatable(df_tab, caption = input$sets_selin, options = list(
    pageLength = 15,
    lengthMenu = list(c(15, 25, 50, 75, 100, -1), c("15", "25", "50", "75", "100", "All"))
  ))
  return(dt)
}

create_central_meas <- function(settings, nms_setting){
  for(i in 1:length(settings)){
    settings[[i]] <- getNZ(settings[[i]], vals$ntraits, input$cor_t, input$cor_m)
    settings[[i]] <- get_subnet(settings[[i]])
    settings[[i]]@x[abs(settings[[i]]@x) > 0] <- 1
  }

  settings <- settings[lapply(settings, function(y) length(y@x)) > 0]

  shiny::validate(shiny::need(length(settings) > 0, "No Connections"))

  graph_mats <- lapply(settings, igraph::graph_from_adjacency_matrix, mode = "undirected")
  graph_mats_connected <- lapply(graph_mats, del_iso_nodes)
  suppressWarnings(graph_dfs <- mapply(function(g, nm) data.frame(igraph::vertex_attr(g, "name"), igraph::degree(g), igraph::betweenness(g), igraph::closeness(g), nm),
                                       graph_mats_connected, nms_setting, SIMPLIFY = FALSE))

  graph_df <- do.call(rbind.data.frame, graph_dfs)

  if(vals$mode == "gxe"){
    sett_label <- "Environment"
  }
  else{
    sett_label <- "Setting"
  }

  colnames(graph_df) <- c("Name", "Degree", "Betweenness", "Closeness", "Setting")

  df_graphs_reshaped <- reshape(data = graph_df,
                                direction = "long",
                                idvar = c("Name", "Setting"),
                                varying = c("Degree", "Betweenness", "Closeness"),
                                times = c("Degree", "Beweenness", "Closeness"),
                                v.names = "Value",
                                timevar = "Centrality")

  if(shiny::isTruthy(vals$map_nodes)){
    df_graphs_reshaped <- merge(df_graphs_reshaped, vals$map_nodes, by.x = "Name", by.y = "node")
    df_graphs_reshaped <- df_graphs_reshaped[order(factor(df_graphs_reshaped$Name, levels = vals$map_nodes$node)),]
    df_graphs_reshaped$Name <- factor(df_graphs_reshaped$Name, levels = vals$map_nodes$node)
    df_graphs_reshaped$node_group <- factor(df_graphs_reshaped$node_group, levels = unique(vals$map_nodes$node_group))

    if(input$cen_meas_col_switch){
      rects <- merge(graph_df, vals$map_nodes, by.x = "Name", by.y = "node")
      rects <- rects[, c(1, 6, 7)]
      rects <- distinct(rects, .keep_all = TRUE)
      rects <- rects[order(factor(rects$Name, levels = vals$map_nodes$node)),]

      a <- tapply(seq_along(rects$node_group), rects$node_group, min)
      a <- a[!is.na(a)]
      a <- data.frame(from = a, node_group = names(a), row.names = NULL)

      b <- tapply(seq_along(rects$node_group), rects$node_group, max)
      b <- b[!is.na(b)]
      b <- data.frame(to = b, node_group = names(b), row.names = NULL)

      rects_com <- merge(a, b, by.x = "node_group", by.y = "node_group")
      rects_com <- merge(rects_com, rects[, 2:3], by.x = "node_group", by.y = "node_group")

      rects_com <- distinct(rects_com, .keep_all = TRUE)
      rects_com <- rects_com[order(factor(rects_com$node_group, levels = unique(vals$map_nodes$node_group))),]
      rects_com$node_group <- factor(rects_com$node_group, levels = unique(vals$map_nodes$node_group))

      rects_com$to <- as.numeric(rects_com$to) + 1
      rects_com$to <- rects$Name[rects_com$to]

      rects_com$from <- rects$Name[rects_com$from]

      rects_com[nrow(rects_com), 3] = rects$Name[nrow(rects)]
    }
  }

  if(input$cen_meas_leg_pos %in% c("bottom", "top")){
    div <- 3
    box <- "horizontal"
  }
  else{
    div <- 8
    box <- "vertical"
  }

  p <-  ggplot2::ggplot() +
    ggplot2::geom_jitter(data = df_graphs_reshaped, ggplot2::aes(x = Name, y = Value, color = Setting), size = 3) +
    ggplot2::facet_grid(Centrality ~ ., scales = "free_y") +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::labs(color = sett_label) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1, size = 13),
                   legend.text = ggplot2::element_text(size = 15),
                   axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title.x = ggplot2::element_text(size = 16, face = "bold"),
                   strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "#F5F2F2"),
                   panel.background = ggplot2::element_rect(fill = "#F5F2F2", color = "black"),
                   legend.background = ggplot2::element_rect(fill = "#F5F2F2"),
                   legend.key = ggplot2::element_rect(fill = "#F5F2F2"),
                   legend.key.size = ggplot2::unit(3, "line"),
                   legend.title =  ggplot2::element_blank(),
                   legend.position = input$cen_meas_leg_pos,
                   legend.box = box)

  if(input$cen_meas_col_switch && shiny::isTruthy(vals$map_nodes)){
    p <- p +
      ggplot2::geom_rect(data = rects_com, ggplot2::aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf, fill = node_group), alpha = 0.2) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = ceiling(length(rects_com$node_group)/div))) +
      ggplot2::scale_fill_manual(values = rects_com$node_color)

  }

  return(p)
}

del_iso_nodes <- function(graph_object){
  isolated <- which(degree(graph_object) == 0)
  graph_object_connected <- igraph::delete.vertices(graph_object, isolated)
  return(graph_object_connected)
}

#Function that gets the plot for the summary statistics plot
#This function takes in no arguments. When called the function returns a
#ggplot object
get_summ_stats_plot <- function(cor_t, cor_m, meas){
  dat <- data.frame(sett = factor(),
                    meas = factor(),
                    val = numeric())
  for (i in vals$sett_names){
    mat <- getNZ(vals$networks[[i]], vals$ntraits, cor_t, cor_m)
    mat <- get_subnet(mat)

    g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
    g2 <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = NULL)
    mat <- del_iso_nodes(g)

    dat <- rbind(dat, data.frame("sett" = i, "meas" ="# of nodes", "val" = igraph::vcount(g)))
    dat <- rbind(dat, data.frame("sett" = i, "meas" = "# of edges", "val" = igraph::ecount(g)))
    dat <- rbind(dat, data.frame("sett" = i, "meas" = "Global Clustering Coefficient", "val" = igraph::transitivity(g)))
    dat <- rbind(dat, data.frame("sett" = i, "meas" = "Average Clustering Coefficient", "val" = igraph::transitivity(g, type = "average")))

    #When users input their own function
    if(meas != ""){
      for(str in unlist(strsplit(meas, ";"))){
        str <- trimws(str, which = "both")
        func <- unlist(strsplit(str, ","))

        if(!exists(func[1], mode = "function")){
          next
        }

        tot_func <- paste0(func[1], "(g")

        if(length(func) > 1){
          for(j in func[2:length(func)]){
            tot_func <- paste0(tot_func, ",", j)
          }
          tot_func <- paste0(tot_func, ")")
        }

        else{
          tot_func <- paste0(tot_func, ")")
        }
        tryCatch(expr = {
          to_eval <- parse(text = tot_func)
          dat <- rbind(dat, data.frame("sett" = i, "meas" = func[[1]], "val" = eval(to_eval)))
        },
        error = function(cond) NULL
        )
        #to_eval <- parse(text = tot_func)
        #try(dat <- rbind(dat, data.frame("sett" = i, "meas" = func[[1]], "val" = eval(to_eval))), silent = TRUE)
      }
    }
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = sett, y = val, fill = sett)) +
    ggplot2::geom_bar(stat = 'identity', width = 0.95) +
    ggplot2::facet_wrap(~meas, ncol = 3, scale = "free_y") +
    ggplot2::labs(fill = NULL) +
    ggplot2::ylab("Value") +
    ggplot2::xlab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                   axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title.x = ggplot2::element_text(size = 16, face = "bold"),
                   strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "#F5F2F2"),
                   panel.background = ggplot2::element_rect(fill = "#F5F2F2"),
                   legend.position = "none")

  return(p)
}

#Helper function for getting the plots for bootstrap.
#This matrix takes in either an adjacency matrix or an partial correlation
#matrix from one bootstrap sample, it also takes a marker, and a boolean
#value (trait) that indicates whether it is a trait or not
#This function then gets all the names of the other markers that the chosen marker
#has a connection with in a particular bootstrap sample
helper_bootstrap_func <- function(mat, mkr){
  lst <- names(which(mat[mkr, ] != 0))
  lst <- lst[lst != mkr]
  return(lst)
}

#This function creates the bootstrap plots
#This function takes no argument and returns the resamples information on the chosen marker
#or trait chosen in the dropdownmenu
bootstrap_func <- function(new_res = NULL){
  mkr <- input$marker
  res_boots <- vector("list", length(futureData$data10))
  for (i in 1:length(futureData$data10)){
    res_boots[[i]] <- unlist(lapply(futureData$data10[[i]], helper_bootstrap_func, mkr))
  }
  data <- stack(setNames(res_boots, seq_along(res_boots)))
  names(data) <- c("node", "environment")

  tab <- table(data$node, data$environment)
  for(i in 1:dim(tab)[2]){
    tab[, i] <- tab[, i]/length(futureData$data10[[i]])
  }
  tab <- data.frame(tab)
  colnames(tab) <- c("node", "environment", "freq")
  tab <- merge(tab, vals$map_nodes, by.x = "node", by.y = "node")
  tab$environment <- factor(tab$environment, levels = as.character(1:length(vals$networks)), labels = vals$sett_names)
  tab <- tab[order(factor(tab$node, levels = vals$map_nodes$node)),]
  tab$node <- factor(tab$node, levels = vals$map_nodes$node)
  tab$node_group <- factor(tab$node_group, levels = unique(vals$map_nodes$node_group))
  row.names(vals$map_nodes) <- vals$map_nodes$node
  frame_colors <- unique(vals$map_nodes[unique(as.character(tab$node)), ]$node_color)
  #We create a dataframe that plots diamonds shapes on top of markers
  #that in the real data also has a connection with marker/trait. True positives.
  df_conn <- data.frame(node = factor(),
                        environment = factor(),
                        freq = as.numeric())

  for (i in vals$sett_names){
    mat <- vals$networks[[i]][, mkr]
    mks_all <- vals$node_names[which(mat != 0)]
    mks_all <- mks_all[mks_all != mkr]
    tab_new <- tab[tab$environment == i, ]
    df_conn <- rbind(df_conn, tab_new[tab_new$node %in% mks_all, 1:3])
  }
  lgnd.pos <- ifelse(is.null(vals$map_nodes), "none", "right")

  if(vals$mode == "gxe"){
    text_tooltip <- "Chromosome: "
  }
  else{
    text_tooltip <- "Group: "
  }
  p <- suppressWarnings(
    ggplot2::ggplot() +
      ggplot2::geom_bar(data = tab, ggplot2::aes(x = node,
                                                 y = freq,
                                                 fill = node_group,
                                                 text = paste0("Name: ",  node, "\n", text_tooltip, node_group),
                                                 text2 = paste0("Percentage", freq * 100)),
                        stat = 'identity') +
      ggplot2::geom_point(data = df_conn, ggplot2::aes(x = node, y = freq), shape = 23, fill = "red", color = "black", color = "red", size = 2) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::facet_wrap(~environment, nrow = length(vals$networks)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1, size = 10),
                     legend.background = ggplot2::element_rect(fill = "#F5F2F2"),
                     legend.text = ggplot2::element_text(size = 15),
                     strip.text = ggplot2::element_text(face = "bold"),
                     axis.title.y = ggplot2::element_text(size = 15),
                     axis.title.x = ggplot2::element_text(size = 16, face = "bold")) +
      ggplot2::labs(fill = NULL) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::scale_fill_manual(values = frame_colors) +
      ggplot2::theme(legend.position = lgnd.pos))

  gp <- plotly::ggplotly(p, tooltip = c("text", "text2"))

  for(i in 1:length(gp$x$data)){
    text <- gp$x$data[[i]]$text
    text <- gsub("paste0(\"Percentage\", freq * 100): Percentage", "Certainty: ", text, fixed = TRUE)
    if(is.null(vals$mapping)){
      if(vals$mode == "gxe"){
        text <- gsub("Chromosome:(.*?)>", "", text)
      }
      else{
        text <- gsub("Group:(.*?)>", "", text)
      }

    }
    gp$x$data[[i]]$text <- text
  }

  gp <- gp %>%
    plotly::layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor = 'rgba(0,0,0,0)')

  return(gp)
}

#Function that takes two matrices as argument and gets the distance between the two matrices
#based on a distance metric. The function can calculate the distance between matrices based
#on the weighted matrices or adjacency matrices
mat_dist <- function(mat1, mat2, dist = "Euclidean", mat_type){
  if(mat_type == "Adjacency"){
    mat1@x[abs(mat1@x) > 0] <- 1
    mat2@x[abs(mat2@x) > 0] <- 1
  }

  diff_mat <- mat1 - mat2

  if(dist == "Euclidean"){
    diff_mat <- diff_mat^2
    all_sum <- sqrt(sum(diff_mat))

  }
  else if(dist == "Manhattan"){
    diff_mat <- abs(diff_mat)
    all_sum <- sum(diff_mat)
  }

  else if(dist == "Canberra"){
    diff_mat <- sum(abs(diff_mat))
    denom <- sum(abs(mat1) + abs(mat2))
    all_sum <- diff_mat/denom
  }

  else if(dist == "Jaccard"){
    intersection <- sum(mat1 == mat2)
    denom <- mat1@dim[[1]]^2 - intersection
    all_sum <- 1 - intersection/denom
  }

  return(all_sum)
}

apply_mat_dist_list <- function(mat_list, dist = "Euclidean", mat_type, for_plot = FALSE){
  names(mat_list) <- vals$sett_names
  combis <- t(combn(names(mat_list), 2))

  if(for_plot){
    combis_df <- data.frame(combis)
    colnames(combis_df) <- c("net1", "net2")
    for(d in dist){
      distances <- apply(combis, 1,
                         function(x) mat_dist(mat_list[[x[1]]], mat_list[[x[2]]], dist = d, mat_type = mat_type))
      combis_df[[paste0(d)]] <- distances
    }

    combis_df <- combis_df[order(factor(combis_df$net1, levels = vals$sett_names)),]
    combis_df$net2 <- factor(combis_df$net2, levels = vals$sett_names)
    combis_df$net1 <- factor(combis_df$net1, levels = vals$sett_names)
    nc <- ncol(combis_df)
    new_combis_df <- combis_df[, c(2, 1, 3:nc)]
    colnames(new_combis_df) <- colnames(combis_df)
    combis_df <- rbind(combis_df, new_combis_df)
    combis_df <- reshape(combis_df,
                         varying = colnames(combis_df[3:nc]),
                         v.names = "value",
                         timevar = "Distance Measure",
                         times = colnames(combis_df[3:nc]),
                         direction = "long")

    if(isTRUE(input$log_check)){
      combis_df$value <- log(combis_df$value)
    }

    p <- ggplot2::ggplot(combis_df, ggplot2::aes(x = factor(net2), y = value, fill = `Distance Measure`)) +
      ggplot2::geom_bar(stat = "identity", width = 0.5, position = "dodge") +
      ggplot2::facet_wrap(~net1, scales = "free_x") +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Value") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
                     axis.title.x = ggplot2::element_text(size = 16, face = "bold"),
                     strip.text = ggplot2::element_text(size = 14, face = "bold"))


    return(p)
  }

  else{
    distances <- apply(combis, 1,
                       function(x) mat_dist(mat_list[[x[1]]], mat_list[[x[2]]], dist = dist, mat_type = mat_type))
    distances_mat <- matrix(nrow = length(mat_list), ncol = length(mat_list))
    distances_mat[lower.tri(distances_mat)] <- distances

    rownames(distances_mat) <- vals$sett_names
    colnames(distances_mat) <- vals$sett_names
    if(ncol(distances_mat) > 1){
      distances_mat <- distances_mat[-1, -ncol(distances_mat)]
    }
  }

  return(distances_mat)
}

comm_detection_plot <- function(setting){
  dat <- data.frame(sett = factor(),
                    meas = factor(),
                    val = numeric())
  set.seed(6708)
  mat <- getNZ(vals$networks[[setting]], vals$ntraits, input$cor_t, input$cor_m)
  mat <- get_subnet(mat)

  shiny::validate(shiny::need(length(mat@x) > 0, "No Connections"))

  mat@x[abs(mat@x) > 0] <- 1
  g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = NULL)
  g <- del_iso_nodes(g)

  lay <- igraph::layout_with_fr(g)
  row.names(lay) <- names(V(g))

  if(input$cluster_algs == "Fast Greedy"){
    try(c1 <- igraph::cluster_fast_greedy(g))
  }
  else if(input$cluster_algs == "Edge Betweenness"){
    try(c1 <- igraph::cluster_edge_betweenness(g))
  }

  for (i in 1:nrow(lay)){
    if(is.element(rownames(lay)[i], rownames(coords()))){
      lay[rownames(lay)[i], ] <- coords()[rownames(lay)[i], ]
    }
  }

  sel_nodes <- dimnames(mat)[[1]]
  if(!is.null(vals$map_nodes)){
    frame_colors <- vals$map_nodes[vals$map_nodes$node %in% sel_nodes, ]$node_color
  }

  p <- plot(c1, g, layout = cbind(lay[, 1], -1 * lay[, 2]), vertex.shape = "fcircle",
            edge.curved = input$roundness)
  return(p)
}

modularity_plot <- function(){
  dat <- data.frame(sett = factor(),
                    meas = factor(),
                    val = numeric())
  for (i in vals$sett_names){
    mat <- getNZ(vals$networks[[i]], vals$ntraits, input$cor_t, input$cor_m)
    mat <- get_subnet(mat)
    mat@x[abs(mat@x) > 0] <- 1
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = NULL)
    g <- del_iso_nodes(g)

    try(dat <- rbind(dat, data.frame("sett" = i, "meas" ="Fast Greedy", "val" = igraph::modularity(igraph::cluster_fast_greedy(g)))))
    try(dat <- rbind(dat, data.frame("sett" = i, "meas" = "Leading Eigenvector", "val" = igraph::modularity(igraph::cluster_leading_eigen(g)))))
    if(input$meas_butt != ""){
      for(str in unlist(strsplit(input$meas_butt, ";"))){
        str <- trimws(str, which = "both")
        func <- unlist(strsplit(str, ","))

        if(!exists(func[1], mode = "function")){
          next
        }

        tot_func <- paste0(func[1], "(g")
        if(length(func) > 1){
          for(j in func[2:length(func)]){
            tot_func <- paste0(tot_func, ",", j)
          }
          tot_func <- paste0(tot_func, ")")
        }

        else{
          tot_func <- paste0(tot_func, ")")
        }
        tot_func <- paste0("modularity(", tot_func, ")")
        to_eval <- parse(text = tot_func)
        dat <- rbind(dat, data.frame("sett" = i, "meas" = func[1], "val" = eval(to_eval)))
      }
    }
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = sett, y = val, fill = sett)) +
    ggplot2::geom_bar(stat = 'identity', width = 0.95) +
    ggplot2::facet_wrap(~meas, ncol = 3, scale = "free_y") +
    ggplot2::ylab("Value") +
    ggplot2::xlab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                   axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title.x = ggplot2::element_text(size = 16, face = "bold"),
                   strip.text = ggplot2::element_text(size = 14, face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "#F5F2F2"),
                   panel.background = ggplot2::element_rect(fill = "#F5F2F2"),
                   legend.position = "none")

  return(p)
}

#Columns that takes in the users inputted data (character vector) for columns to exclude and returns
#a numeric vector with the columns to remove in it
#This function takes in an optional argument (whether to return warning or not) and returns
#a numeric vector
col_to_exclude <- function(vec, warn = FALSE){
  suppressWarnings(vec_col <- as.numeric(trimws(unlist(strsplit(vec," ")))))
  if(length(vec_col) > 0){
    if(sum(is.na(vec_col)) > 0){
      if(isTRUE(warn)) shiny::showNotification("Input given is not only numbers", type = "warning")
      vec_col <- vec_col[!is.na(vec_col)]
    }
    if(sum(vec_col < 1) > 0){
      if(isTRUE(warn)) shiny::showNotification("A number in the input is less than 1", type = "warning")
      vec_col <- vec_col[vec_col > 0]
    }
    if(!all(vec_col == floor(vec_col))){
      if(isTRUE(warn)) shiny::showNotification("A non integer was given", type = "warning")
      vec_col <- vec_col[vec_col == floor(vec_col)]
    }
  }
  return(vec_col)
}

#This function checks if the column numbers to exclude are larger than the amount of columns
#This function takes in a dataframe and a numeric vector as arguments
#and returns a dataframe with potentially less columns
check_cols_to_excl <- function(df, vec){
  big_vec <- vec > ncol(df)
  if(shiny::isTruthy(big_vec) && sum(big_vec) > 0){
    shiny::showNotification(paste("Some numbers passed are bigger than the amount of columns in data",
                                  "(", toString(vec[big_vec]), ")", sep = ""), type = "warning")
    vec <- vec[!big_vec]
  }
  if(shiny::isTruthy(vec) && length(vec) > 0){
    if((ncol(df) - length(vec)) < 2){
      shiny::showNotification("Dataframe will be left with less than two columns", type = "error")
      shiny::updateNumericInput(session = session, inputId = "exc_columns_mapping", value = "")
      return(df)
    }
    else{
      return(df[, -vec])
    }
  }
  else{
    return(df)
  }
}

#Function that performs networks reconstruction in the startup modals
#This function takes in no arguments and returns a list with the partial correlation
#matrices with the reconstructed networks from the uploaded dataframes
perform_startup_recon <- function(l_args) {
  len_files <- length(uploadedFiles$files)
  nms <- names(uploadedFiles$files)
  if(is.null(vals$networks)){
    vals$networks <- list()
  }
  shiny::withProgress(message = "Reconstructing Networks", value = 0, {
    for (i in 1:len_files) {
      curr_name <- tools::file_path_sans_ext(nms[[i]])
      shiny::incProgress(1/len_files, detail = paste("Reconstructing Network", curr_name))
      if(curr_name %in% names(vals$networks)){
        shiny::showNotification(paste("File with name ", curr_name, "has already been reconstructed. \n", "Skipping."), type = "warning")
        next
      }
      tryCatch(
        expr = {
          vals$networks[[curr_name]] <- netgwas::netphenogeno(
            data = uploadedFiles$files[[i]],
            method = l_args[["net_method_start"]],
            rho = l_args[["net_rho_start"]],
            n.rho = l_args[["net_n.rho_start"]],
            rho.ratio = l_args[["net_rho.ratio_start"]],
            ncores = l_args[["net_ncores_start"]],
            em.iter = l_args[["net_em.iter_start"]],
            em.tol = l_args[["net_em.tol_start"]],
            verbose = FALSE
          )
          vals$networks[[curr_name]] <- netgwas::selectnet(
            netgwas.obj = vals$networks[[curr_name]],
            opt.index = l_args[["sel_opt.index_start"]],
            criteria = l_args[["sel_criteria_start"]],
            ebic.gamma = l_args[["sel_ebic.gamma_start"]],
            ncores = l_args[["sel_ncores_start"]],
            verbose = FALSE
          )$par.cor
        },
        error = function(cond) {
          shiny::showNotification(
            paste("File ", nms[[i]], " was unsuccesful during reconstruction"),
            type = "error",
            duration = NULL
          )
        }
      )
    }
  })
}

#Function that perform some processing on the netphengeno and selectnet arguments if needed
#This function takes in the argument that needs processing and a boolean that says
#whether the argument is for ncores or not
process_args_recon <- function(arg, nc){
  if(isTRUE(nc)){
    if(arg == "all"){
      return(all)
    }
    else{
      x <- tryCatch(expr = {as.numeric(arg)},
                    warning = function(cond) return(1),
                    error = function(cond) return(1))
      return(x)
    }
  }
  else{
    if(arg == "NULL"){
      return(NULL)
    }
    else{
      x <- tryCatch(expr = {as.numeric(arg)},
                    warning = function(cond) return(NULL),
                    error = function(cond) return(NULL))
      return(x)
    }
  }
}

#Function to get the arguments for the netphenogeno and selectnet function in the proper format
#This function takes a boolean argument and returns a list with all of the arguments for the functions
get_args_recon <- function(start_up){
  args_net <- args_netphenogeno[-length(args_netphenogeno)]
  args_sel <- args_selectnet[-length(args_selectnet)]

  if(isTRUE(start_up)){
    args_net_mod <- paste("net_", args_net, "_start", sep = "")
    args_sel_mod <- paste("sel_", args_sel, "_start", sep = "")
    args_all_mod <- c(args_net_mod, args_sel_mod)
    check1 <- c("net_rho_start", "sel_opt.index_start", "sel_criteria_start")
    check2 <- c("net_ncores_start", "sel_ncores_start")
  }
  else{
    args_net_mod <- paste("net_", args_net, sep = "")
    args_sel_mod <- paste("sel_", args_sel, sep = "")
    args_all_mod <- c(args_net_mod, args_sel_mod)
    check1 <- c("net_rho", "sel_opt.index", "sel_criteria")
    check2 <- c("net_ncores", "sel_ncores")
  }

  #NEED TO USE LAPPLY INSTEAD OF FOR LOOP

  list_args <- setNames(vector("list", length = 11), args_all_mod)
  for(i in 1:length(args_all_mod)){
    if(args_all_mod[i] %in% check1){
      list_args[[args_all_mod[i]]] <- lapply(input[[args_all_mod[[i]]]], process_args_recon, FALSE)[[1]]
    }
    else if(args_all_mod[i] %in% check2){
      list_args[[args_all_mod[i]]] <- lapply(input[[args_all_mod[[i]]]], process_args_recon, TRUE)[[1]]
    }
    else{
      list_args[[args_all_mod[i]]] <- input[[args_all_mod[i]]]
    }
  }
  return(list_args)
}


map_nodes_to_group <- function(trt_typs){
  if(isTRUE(input$gxe_mode) && !shiny::isTruthy(vals$map_nodes)){
    n_markers <- length(vals$node_names) - vals$n_traits
    vals$map_nodes <- data.frame("node" = vals$node_names, "node_group" = c(rep("Trait", vals$n_traits), rep("Marker", n_markers)))
  }

  if(isTruthy(input$trait_types)){
    trt_grouping <- data.frame("node" = vals$node_names[1:vals$n_traits], "node_group"= rep(trt_typs$envs, trt_typs$freq))
    vals$map_nodes$node <- as.character(vals$map_nodes$node)
    vals$map_nodes$node_group <- as.character(vals$map_nodes$node_group)
    if(!sum(trt_grouping$node %in% vals$map_nodes$node)){
      vals$map_nodes <- rbind(trt_grouping, vals$map_nodes)
    }
    else{
      vals$map_nodes[match(trt_grouping$node, vals$map_nodes$node, nomatch = 1), ][, 1:2] <- trt_grouping
    }
  }
}

complete_df <- function(df){
  #If less unique values in first column, then first column contains groups for nodes, or vice versa
  df <- as.data.frame(df)
  node_group <- df$node_group[match(vals$node_names, df$node)]
  node_group <- as.character(node_group)
  node_group[is.na(node_group)] <- ifelse(vals$mode == "gxe", "Trait", "na")
  new_df <- data.frame("node" = vals$node_names, "node_group" = factor(node_group))
  node_color <- custom.col[1:length(levels(new_df$node_group))]
  node_color <- node_color[match(new_df$node_group, levels(new_df$node_group))]
  new_df$node_color <- node_color
  return(new_df)
}
