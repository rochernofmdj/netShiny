# Not in operator
`%nin%` <- Negate(`%in%`)

# Custom colors that we (semi) manually selected to use in the plots
# These colors were chosen such that any two colors next to each are as distinct as possible
custom.col <- c("#DBB165",  "#52854C", "#4E84C4", "#C3D7A4","#C4961A", "#8B4513", "#293352",
                "#E69F00", "#56B4E9", "#009E73", "#B62A3D", "#0072B2", "#D55E00", "#CC79A7",
                "#000000", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
                "#B15928", "#FFD700", "#C5B358")

###WARNING: #009E73 AND #1B9E77 ARE REALLY SIMILAR TO EACH OTHER

# Custom bsModal that does not close when users click outside the dialog modal
# or when users press the escape key
custombsModal <- function(...) {
  bmodal <- shinyBS::bsModal(...)
  bmodal[[2]]$`data-backdrop` <- "static"
  bmodal[[2]]$`data-keyboard` <- "false"
  return(bmodal)
}

# This function checks if only a few nodes are chosen to create a subnetwork
# If a subnetwork needs to be created, it does so by only getting the nodes that are chosen
# with their neighbors
get_subnet <- function(vals, mat){
  if(!is.null(vals$subgraph_nodes)){
    ind_subgraph_nodes <- which(dimnames(mat)[[2]] %in% vals$subgraph_nodes)
    summ_df <- Matrix::summary(mat)
    summ_df <- summ_df[summ_df$i %in% ind_subgraph_nodes | summ_df$j %in% ind_subgraph_nodes, ]
    nms_subgraph_nodes <- sort(union(summ_df$i, summ_df$j))
    nms_subgraph_nodes <- dimnames(mat)[[2]][nms_subgraph_nodes]
    mat <- mat[nms_subgraph_nodes, nms_subgraph_nodes]
  }
  return(mat)
}

# Function that returns an adjacency or weighted matrix with only nodes
# that satisfy threshold value(s) given. If (both) threshold value(s) are NULL, it will
# return an adjacency or weight matrix that only include non-isolated nodes
getNZ <- function(vals, input, mat, diff = FALSE, to_plot = FALSE){
  shiny::validate(shiny::need(!is.null(mat), "Getting proper data"))
  if(vals$mode == "gxe"){
    shiny::validate(shiny::need(!is.null(vals$n_traits), "Getting proper data"))
  }
  mat <- Matrix::Matrix(mat)
  n_nodes <- dim(mat)[[2]]
  if(isFALSE(vals$mode == "gxe") || isTRUE(diff) || isTRUE(to_plot)){
    if(!is.null(input$cor_t) && isFALSE(diff)){
      #mat[abs(mat) < input$cor_t] <- 0
      indic <- Matrix::which(abs(mat) < input$cor_t, arr.ind = TRUE)
      mat[indic] <- 0
    }
    diag(mat) <- 0 #Makes diagonal values 0
    mat <- Matrix::drop0(mat) #Returns sparse matrix with no explicit zeroes (including removing diagonal zeroes)
    nzvec <- which(Matrix::colSums(mat) == 0) #Gets indices of isolated nodes
    if(length(nzvec) > 0) mat <- mat[-nzvec, -nzvec]
    invisible(mat)
  }

  else{
    if (!is.null(input$cor_m)){
      mat[(vals$n_traits+1):n_nodes, (vals$n_traits+1):n_nodes][abs(mat[(vals$n_traits+1):n_nodes, (vals$n_traits+1):n_nodes]) < input$cor_m] <- 0
    }
    if (!is.null(input$cor_t)){
      mat[, 1:vals$n_traits][abs(mat[, 1:vals$n_traits]) < input$cor_t] <- 0
      mat[1:vals$n_traits, ][abs(mat[1:vals$n_traits, ]) < input$cor_t] <- 0
    }
    diag(mat) <- 0  #Makes diagonal values 0
    mat <- Matrix::drop0(mat, tol = 0) #Returns sparse matrix with no explicit zeroes for opt.adj, thus removing diagonal 0 values
    nzvec <- which(Matrix::colSums(mat) == 0) #gets indices of traits and markers that do not have any links
    nzvec <- nzvec[nzvec > vals$n_traits] #gets indices of only markers that do not have any links
    mat <- mat[-nzvec, -nzvec]
    invisible(mat)
  }
}

# This function takes as argument an adjacency or weight matrix
# The function returns the names of only non-isolated nodes
get_nz_nodes <- function(mat, vals, input){
  mat <- getNZ(vals = vals, input = input, mat = mat)
  nz_names <- dimnames(mat)[[2]]
  return(nz_names)
}

del_iso_nodes <- function(graph_object){
  isolated <- which(igraph::degree(graph_object) == 0)
  graph_object_connected <- igraph::delete.vertices(graph_object, isolated)
  return(graph_object_connected)
}
