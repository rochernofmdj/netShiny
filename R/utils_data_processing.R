# Arguments of netphenogeno and selectnet
args_netphenogeno <- methods::formalArgs(netgwas::netphenogeno)[-1]
args_netphenogeno <- args_netphenogeno[-length(args_netphenogeno)]
args_selectnet <- methods::formalArgs(netgwas::selectnet)[-1]
args_selectnet <- args_selectnet[-length(args_selectnet)]

# Columns that takes in the users inputted data (character vector) for columns to exclude and returns
# a numeric vector with the columns to remove in it
# This function takes in an optional argument (whether to return warning or not) and returns
# a numeric vector
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

# This function checks if the column numbers to exclude are larger than the amount of columns
# This function takes in a dataframe and a numeric vector as arguments
# and returns a dataframe with potentially less columns
check_cols_to_excl <- function(session, df, vec){
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

# Function that performs networks reconstruction in the startup modals
# This function takes in no arguments and returns a list with the partial correlation
# matrices with the reconstructed networks from the uploaded dataframes
perform_startup_recon <- function(val_nets, files, l_args) {
  len_files <- length(files)
  nms <- names(files)
  if(is.null(val_nets)){
    val_nets <- list()
  }
  shiny::withProgress(message = "Reconstructing Networks", value = 0, {
    for (i in 1:len_files) {
      curr_name <- paste0(tools::file_path_sans_ext(nms[i]), "_", l_args[["net_method_start"]])
      shiny::incProgress(1/len_files, detail = paste("Reconstructing Network", curr_name))
      if(curr_name %in% names(val_nets)){
        shiny::showNotification(paste("File with name ", curr_name, "has already been reconstructed. \n", "Skipping."), type = "warning")
        next
      }
      tryCatch(
        expr = {
          val_nets[[curr_name]] <- netgwas::netphenogeno(
            data = files[[i]],
            method = l_args[[nms[i]]][["net_method_start"]],
            rho = l_args[[nms[i]]][["net_rho_start"]],
            n.rho = l_args[[nms[i]]][["net_n.rho_start"]],
            rho.ratio = l_args[[nms[i]]][["net_rho.ratio_start"]],
            ncores = l_args[[nms[i]]][["net_ncores_start"]],
            em.iter = l_args[[nms[i]]][["net_em.iter_start"]],
            em.tol = l_args[[nms[i]]][["net_em.tol_start"]],
            verbose = FALSE
          )
          val_nets[[curr_name]] <- netgwas::selectnet(
            netgwas.obj = val_nets[[curr_name]],
            opt.index = l_args[[nms[i]]][["sel_opt.index_start"]],
            criteria = l_args[[nms[i]]][["sel_criteria_start"]],
            ebic.gamma = l_args[[nms[i]]][["sel_ebic.gamma_start"]],
            ncores = l_args[[nms[i]]][["sel_ncores_start"]],
            verbose = FALSE
          )$par.cor
        },
        error = function(cond) {
          shiny::showNotification(
            paste("File ", nms[i], " was unsuccesful during reconstruction"),
            type = "error",
            duration = NULL
          )
        }
      )
    }
  })

  return(val_nets)
}

# Function that perform some processing on the netphengeno and selectnet arguments if needed
# This function takes in the argument that needs processing and a boolean that says
# whether the argument is for ncores or not
process_args_recon <- function(arg, nc){
  if(isTRUE(nc)){
    if(arg == "all"){
      return("all")
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

# Function to get the arguments for the netphenogeno and selectnet function in the proper format
# This function takes a boolean argument and returns a list with all of the arguments for the functions
get_args_recon <- function(input, start_up){
  if(isTRUE(start_up)){
    args_net_mod <- paste("net_", args_netphenogeno, "_start", sep = "")
    args_sel_mod <- paste("sel_", args_selectnet, "_start", sep = "")
    args_all_mod <- c(args_net_mod, args_sel_mod)
    check1 <- c("net_rho_start", "sel_opt.index_start")
    check2 <- c("net_ncores_start", "sel_ncores_start")
  }
  else{
    args_net_mod <- paste("net_", args_netphenogeno, sep = "")
    args_sel_mod <- paste("sel_", args_selectnet, sep = "")
    args_all_mod <- c(args_net_mod, args_sel_mod)
    check1 <- c("net_rho", "sel_opt.index", "sel_criteria")
    check2 <- c("net_ncores", "sel_ncores")
  }

  list_args <- stats::setNames(vector("list", length = 11), args_all_mod)
  for(i in 1:length(args_all_mod)){
    curr_arg <- args_all_mod[i]
    if(curr_arg %in% check1){
      list_args[[args_all_mod[i]]] <- process_args_recon(input[[curr_arg]], FALSE)
    }
    else if(curr_arg %in% check2){
      list_args[[curr_arg]] <- process_args_recon(input[[curr_arg]], TRUE)
    }
    else{
      list_args[[curr_arg]] <- input[[curr_arg]]
    }
  }

  return(list_args)
}


map_nodes_to_group <- function(vals, input, trt_typs) {
  map_nodes <- vals$map_nodes
  if (isTRUE(input$gxe_mode) && !shiny::isTruthy(map_nodes)) {
    markers <- setdiff(vals$node_names, trt_typs$node)
    map_nodes <- rbind(trt_typs, data.frame("node" = markers, "node_group" = rep("Marker", length(markers))))
    #map_nodes <- data.frame("node" = vals$node_names, "node_group" = c(rep("Trait", vals$n_traits), rep("Marker", n_markers)))
  }

  else if (isFALSE(input$gxe_mode) && !shiny::isTruthy(map_nodes)) {
    no_groups <- setdiff(vals$node_names, trt_typs$node)
    map_nodes <- rbind(trt_typs, data.frame("node" = no_groups, "node_group" = rep("none", length(no_groups))))
  }

  else if(isTruthy(map_nodes)) {
    map_nodes <- map_nodes[, 1:2]
    bool_group_col <- length(unique(map_nodes[, 1])) < length(unique(map_nodes[, 2]))
    if (bool_group_col) {
      colnames(map_nodes)[1:2] <- c("node_group", "node")
    } else {
      colnames(map_nodes)[1:2] <- c("node", "node_group")
    }

    if (isTRUE(input$gxe_mode)) {
      markers <- setdiff(map_nodes$node, trt_typs$node)
      map_nodes <- rbind(trt_typs, data.frame("node" = markers, "node_group" = map_nodes$node_group[match(markers, map_nodes$node)]))
    } else {
      no_groups <- setdiff(vals$node_names, trt_typs$node)
      map_nodes <- rbind(trt_typs, data.frame("node" = no_groups, "node_group" = map_nodes$node_group[match(no_groups, map_nodes$node)]))
    }
    miss_nodes <- setdiff(vals$node_names, map_nodes$node)
    if (length(miss_nodes) > 0) {
      map_nodes <- rbind(map_nodes, data.frame("node" = miss_nodes, "node_group" = rep("na", length(miss_nodes))))
    }
  }

  # if(isTruthy(input$trait_types)){
  #   trt_grouping <- data.frame("node" = vals$node_names[1:vals$n_traits], "node_group"= rep(trt_typs$envs, trt_typs$freq))
  #   map_nodes$node <- as.character(map_nodes$node)
  #   map_nodes$node_group <- as.character(map_nodes$node_group)
  #   if(!sum(trt_grouping$node %in% map_nodes$node)){
  #     map_nodes <- rbind(trt_grouping, map_nodes)
  #   }
  #   else{
  #     map_nodes[match(trt_grouping$node, map_nodes$node, nomatch = 1), ][, 1:2] <- trt_grouping
  #   }
  # }
  return(map_nodes)
}

complete_df <- function(vals){
  # If less unique values in first column, then first column contains groups for nodes, or vice versa
  #df <- as.data.frame(df)
  new_df <- as.data.frame(vals$map_nodes)
  new_df$node_group <- factor(new_df$node_group)
  node_color <- custom.col[1:length(levels(new_df$node_group))]
  node_color <- node_color[match(new_df$node_group, levels(new_df$node_group))]
  new_df$node_color <- node_color

  # node_group <- df$node_group[match(vals$node_names, df$node)]
  # node_group <- as.character(node_group)
  # node_group[is.na(node_group)] <- "na"
  # new_df <- data.frame("node" = vals$node_names, "node_group" = factor(node_group))
  # node_color <- custom.col[1:length(levels(new_df$node_group))]
  # node_color <- node_color[match(new_df$node_group, levels(new_df$node_group))]
  # new_df$node_color <- node_color
  return(new_df)
}

# Function that gets the trait types grouping
# This function takes in a vector and returns a dataframe
# with a trait grouping and its respective total number
get_trait_groups <- function(trait_types){
  # Split the trait types first by comma
  # Then split them by semicolon, and convert into a dataframe, then back into list
  # So we get a list with the grouping with their frequency as a list
  trt_typs <- trimws(strsplit(trait_types, split = ",")[[1]])
  trt_typs_bool <- gsub("\\s", "", trt_typs)
  bools <- grepl("^\\w+:{1}\\d+$", trt_typs_bool)
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

# utils::globalVariables(c("net2", "value", "Distance Measure", "node", "freq", "node_group", "%>%",
#                          "V", "Name", "Value", "Setting", "from", "to", "degree", "weights", "type", "stat",
#                          "width", "density", "sett", "val", "cons", "isTruthy", "HTML", "textOutput",
#                          "myFuture", '<<-', "head", "as_data_frame"))
