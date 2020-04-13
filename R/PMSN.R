#' PM_SN: package for extraction and analysis resource social graph
#'
#' @name pmsn
#' @docType package

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Save social graph into csv file
#' @param gr the igraph object which will be saved in csv file
#' @param res_graph_filename the resulting filename
#' @export
sn_to_csv <- function(gr, res_graph_filename) {
  v_e_list <- igraph::as_data_frame(gr, what = "both")
  df <- v_e_list$edges
  df$width <- df$freq * 0.01
  df$avg_time <- round((as.numeric(df$total_time) / df$freq), 4)
  df$max_time <- round(as.numeric(df$max_time), 4)
  df$min_time <- round(as.numeric(df$min_time), 4)
  df$total_time <- round(as.numeric(df$total_time), 4)
  v_list <- v_e_list$vertices
  df$f_group <- 0
  df$f_leader <- 0.0
  df$f_system_group <- 0
  df$t_group <- 0
  df$t_leader <- 0.0
  df$t_system_group <- 0
  for (i in c(1:nrow(v_list))) {
    if (nrow(df[df$from == v_list$name[i], ]) > 0) {
      if ("group" %in% colnames(v_list)) {
        df[df$from == v_list$name[i], ]$f_group <- v_list$group[i]
        df[df$to == v_list$name[i], ]$t_group <- v_list$group[i]
        }
      if ("leader" %in% colnames(v_list)) {
        df[df$from == v_list$name[i], ]$f_leader <- v_list$leader[i]
        df[df$to == v_list$name[i], ]$t_leader <- v_list$leader[i]
        }
      if ("system_group" %in% colnames(v_list)) {
        df[df$from == v_list$name[i], ]$f_system_group <- v_list$system_group[i]
        df[df$to == v_list$name[i], ]$t_system_group <- v_list$system_group[i]
      }
    }
  }
  if ("leader" %in% colnames(df)) {
    df$f_size <- df$f_leader * 1000
    df[df$f_size < 10, ]$f_size <- 10
    df$t_size <- df$t_leader * 1000
    df[df$t_size < 10, ]$t_size <- 10
  }
  df <- tidyverse::rename(df, source = from, target = to)
  write.csv(df, file = res_graph_filename, row.names = FALSE, quote = FALSE)
}

#' Sends social graph to interactive Cytoscape panel
#' @param gr the igraph object which will be send to Cytoscape interactive panel
#' @export
sn_in_cytpscape <- function(gr) {
  gr <- igraph::set_edge_attr(gr, name = "width", index = E(gr), value = (E(gr)$freq * 0.01))
  E(gr)$avg_time <- round((as.numeric(E(gr)$total_time) / E(gr)$freq), 4)
  E(gr)$max_time <- round(as.numeric(E(gr)$max_time), 4)
  E(gr)$min_time <- round(as.numeric(E(gr)$min_time), 4)
  E(gr)$total_time <- round(as.numeric(E(gr)$total_time), 4)

  gr <- igraph::set_vertex_attr(gr, name = "id", index = V(gr), value = c(1:vcount(gr)))
  size <- c(1:(igraph::vcount(gr)))
  size <- size * 0 + 10
  if ("leader" %in% vertex_attr_names(gr)) {
    size <- V(gr)$leader * 1000
  }
  gr <- igraph::set_vertex_attr(gr, name = "size", index = V(gr), value = size)
  V(gr)$size[V(gr)$size < 10] <- 10

  RCy3::createNetworkFromIgraph(gr,  title = "SN_PM social graph", collection = "SN_PM graphs")

  #Set properties
  RCy3::setVisualStyle("BioPAX_SIF")
  RCy3::layoutNetwork("force-directed-cl")

  #Nodes
  RCy3::setNodeColorDefault("#C0C0C0")
  RCy3::setNodeFillOpacityBypass(V(gr)$name, 225)
  RCy3::setNodeFontSizeDefault(16)
  RCy3::setNodeSizeBypass(V(gr)$name, V(gr)$size)

  #Edges
  e_list <- igraph::as_data_frame(gr, what = "edge")
  RCy3::setEdgeTargetArrowShapeDefault("ARROW")
  RCy3::setEdgeLineWidthBypass(paste(e_list$from, "(interacts with)", e_list$to), e_list$width)
  RCy3::setEdgeLabelBypass(paste(e_list$from, "(interacts with)", e_list$to), e_list$avg_time)
  RCy3::setEdgeFontSizeDefault(14)
}

#' Add system group into social graph as attribute
#' @param gr the igraph object
#' @param nodes_desc_fn filename with mapping nodes and system groups
#' @return igraph object with new attribute
#' @export
add_system_group <- function(gr, nodes_desc_fn) {
  nodes_desc <- read.csv2(file = nodes_desc_fn, header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE, encoding = "en_US.UTF-8")
  system_gr <- integer(vcount(gr))
  for (i in c(1:vcount(gr))) {
    t <- min(nodes_desc[nodes_desc$name == V(gr)$name[i], ]$membership)
    if (length(t) > 0 & !is.na(t) & !is.infinite(t)) {
      system_gr[i] <- t
      }
  }
  gr <- igraph::set_vertex_attr(gr, name = "system_group", index = V(gr), value = system_gr)
  return(gr)
}

#' Filters event log
#' @param event_log data.frame with event log
#' @param mode string from set ("start", "end", "resource", "activity", "case")
#' @param value value in string, it is used for filtering
#' @return data.frame with event log
#' @export
event_log_filtering <- function(event_log, mode, value) {
  if (mode == "start") {
    event_log <- event_log[event_log$timestamp >= as.POSIXct(value, tz = "", format = "%Y-%m-%d %H:%M:%OS"), ]
  }else if (mode == "end") {
    event_log <- event_log[event_log$timestamp <= as.POSIXct(value, tz = "", format = "%Y-%m-%d %H:%M:%OS"), ]
  }else if (mode == "resource") {
    event_log <- event_log[event_log$resource != value, ]
  }else if (mode == "activity") {
    event_log <- event_log[event_log$activity != value, ]
  }else if (mode == "case") {
    event_log <- event_log[event_log$case_id != value, ]
  }
  return(event_log)
}

#' Show event log statistics
#' @param log data.frame with event log
log_info <- function(log) {
  print(paste("In log ", nrow(log), " records"))
  print(paste("In log ", length(unique(log$case_id)), " cases"))
  print(paste("In log ", length(unique(log$activity)), " different activities"))
  print(paste("In log ", length(unique(log$resource)), " different resources"))
  print(paste("Log timeline: ", min(log$timestamp), " - ", max(log$timestamp)))
}

#' Import event log from file
#' File format: csv file with columns: case_id, activity, timestamp, resource
#' @param event_log_path file with event log
#' @return data.frame with event log
#' @export
import_event_log <- function(event_log_path) {
  log <- read.csv2(file = event_log_path, header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE, encoding = "en_US.UTF-8")
  if (!("case_id" %in% colnames(log))) {
    print("Event log must contain case_id column!")
    return()
  }else if (!("activity" %in% colnames(log))) {
    print("Event log must contain activity column!")
    return()
  }else if (!("timestamp" %in% colnames(log))) {
    print("Event log must contain timestamp column!")
    return()
  }else if (!("resource" %in% colnames(log))) {
    print("Event log must contain resource column!")
    return()
  }

  log$resource <- gsub(",", "", log$resource)
  log$timestamp <- gsub("T", " ", log$timestamp)
  log$timestamp <- gsub("Z", "", log$timestamp)
  log$timestamp <- as.POSIXct(log$timestamp, tz = "", format = "%Y-%m-%d %H:%M:%OS")
  log_info(log)
  return(log)
}

#' Built resource social graph as igraph object
#' @param log data.frame with event log
#' @param filename_sn_list temporary filename, for keeping whole social graph (not optimized)
#' @return igraph object with resource social graph
#' @export
social_graph_from_event_log <- function(log, filename_sn_list) {
  temp_log <- data.frame(log$case_id, log$timestamp, log$resource, stringsAsFactors = FALSE)
  temp_log <- unique(temp_log)
  names(temp_log) <- c("case_id", "timestamp", "resource")
  entity_list <- unique(temp_log$case_id)
  res <- NULL
  res$resource1[1] <- "0"
  res$resource2[1] <- "0"
  res$freq[1] <- 0
  res$total_time[1] <- 0
  res$max_time[1] <- 0
  res$min_time[1] <- 0
  res <- data.frame(res, stringsAsFactors = FALSE)
  for (i in entity_list) {
    trace <- temp_log[temp_log$case_id == i, ]
    trace <- trace[order(trace$timestamp), ]
    pred_resource <- trace$resource[1]
    pred_time <- trace$timestamp[1]
    for (j in c(2:nrow(trace))) {
      duration <- as.numeric(abs(difftime(trace$timestamp[j], pred_time, units = "mins")))
      exists_ind <- which(res$resource1 == pred_resource & res$resource2 == trace$resource[j])
      if (length(exists_ind) > 0) {
        res$freq[exists_ind] <- res$freq[exists_ind] + 1
        res$total_time[exists_ind] <- res$total_time[exists_ind] + duration
        res$max_time[exists_ind] <- max(res$max_time[exists_ind], duration)
        res$min_time[exists_ind] <- min(res$min_time[exists_ind], duration)
      }else{
        res <- rbind(res, list(pred_resource, trace$resource[j], 1, duration, duration, duration))
      }
      pred_resource <- trace$resource[j]
      pred_time <- trace$timestamp[j]
    }
  }
  res$total_time <- round(res$total_time, 2)
  res$max_time <- round(res$max_time, 2)
  res$min_time <- round(res$min_time, 2)
  write.csv(res, file = filename_sn_list, row.names = FALSE, quote = TRUE)
  print(paste("social_net ready, find it in a file", filename_sn_list))
  gr <- igraph::graph_from_data_frame(res, directed = TRUE)
  print(gr)
  return(gr)
}

#' Find level for optimising size of resource social graph
#' @param links data.frame with list of graph's edges with thier weight
#' @return numeric value of level
get_level <- function(links) {
  return(round(mean(links$freq)))
}

#' Optimize resource social graph removing rare edges
#' @param sn_path temporary filename, for keeping whole social graph (not optimized)
#' @param level the minimum frequency for optimized resource socail graph edges (optional)
#' @param selfloop boolean parameter show whatever social graph should contain selfloops or not (optional)
#' @return igraph object with resource social graph
#' @export
optimize_graph <- function(sn_path, level = NA, selfloop = TRUE) {
  links <- read.csv2(file = sn_path, header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE, encoding = "en_US.UTF-8")
  if (is.na(level)) {
    level <- get_level(links)
  }
  links <- links[links$freq >= level, ]
  if (!selfloop) {
    links <- links[links$resource1 != links$resource2, ]
  }
  gr <- igraph::graph_from_data_frame(links, directed = TRUE)
  print(gr)
  return(gr)
}

#' Find threshold for internal community links' power
#' @param dist_matrix distance matrix for social graph's nodes, in cells distance metric
#' @return numeric value of threshold
get_community_threshold <- function(dist_matrix) {
  uvec <- sort(unique(as.vector(dist_matrix)))
  hist <- 0
  for (i in c(1:length(uvec))) {
    hist[i] <- sum(dist_matrix == uvec[i])
  }
  return(uvec[which(hist == max(hist[2:(length(hist) - 1)]))])
}

#' Detects social communities in resource social graph
#' @param gr igraph object with resource social graph
#' @param threshold numeric value of threshold for internal community links' power (optional)
#' @return igraph object with resource social graph, new node atrribute for community label
#' @export
edgebased_community_detection <- function(gr, threshold = NA) {
  v_cnt <- vcount(gr)
  dist_matrix <- matrix((c(1:(v_cnt ^ 2)) * 0), nrow = v_cnt, ncol = v_cnt)
  neihgbours <- sapply(V(gr), FUN = function(v) ego(gr, order = 1, nodes = v, mode = c("out"), mindist = 1))
  for (i in c(1:(v_cnt - 1))) {
    for (j in c((i + 1):v_cnt)) {
      dist <- (1 - ((length(intersect(neihgbours[[i]], neihgbours[[j]]))) / (max(length(neihgbours[[i]]), length(neihgbours[[j]])))))
      if (is.na(dist)) {
        dist <- 1
      }
      dist_matrix[i, j] <- dist
      dist_matrix[j, i] <- dist
    }
  }
  if (is.na(threshold)) {
    par(mfrow = c(2, 1))
    threshold <- get_community_threshold(dist_matrix)
    print(paste("Auto threshold:", threshold))
  }
  hist(dist_matrix, breaks = 20, freq = TRUE, xlab = "Distance between nodes", ylab = "Frequency", main = "", cex.lab = 1.5, cex.axis = 1.5)
  par(mfrow = c(1, 1))
  adj_matrix <- dist_matrix
  adj_matrix[adj_matrix > threshold] <- 1
  adj_matrix <- 1 - adj_matrix
  comm_gr <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = "weight", diag = FALSE)
  gr_components <- igraph::components(comm_gr)
  gr <- igraph::set_vertex_attr(gr, name = "group", index = V(gr), value = gr_components$membership)
  print(paste("communities count", gr_components$no))
  print(paste("min power of community", min(gr_components$csize)))
  print(paste("median power of community", median(gr_components$csize)))
  print(paste("avg power of community", mean(gr_components$csize)))
  print(paste("max power of community", max(gr_components$csize)))
  return(gr)
}


#' Identifies leader in resource social graph
#' @param gr igraph object with resource social graph
#' @return igraph object with resource social graph, new node atrribute for leader coefficient
#' @export
find_leader <- function(gr) {
  pg <- igraph::page_rank(gr)
  gr <- igraph::set_vertex_attr(gr, name = "leader", index = V(gr), value = pg$vector)
  pg_df <- data.frame(V(gr)$name, V(gr)$leader)
  pg_df <- pg_df[order(pg_df$V.gr..leader, decreasing = TRUE), ]
  print(head(pg_df))
  return(gr)
}

#' @example
#' run<-function() {
#'   input_filename<-"res_event_log_4500.csv"
#'   sn_list_filename<-"social_net_edge_list.csv"
#'   sn_ve_list_filename<-"social_net_edge_vert_list.csv"
#'   zendesk_group_filename <- "user_desc.csv"
#'   log<-import_event_log(input_filename)
#'   gr<-social_graph_from_event_log(log,sn_list_filename)
#'   summary(E(gr)$freq)
#'   hist(E(gr)$freq)
#'   gr<-optimize_graph(sn_list_filename, selfloop=FALSE)
#'   gr<-find_leader(gr)
#'   e_gr<-edges_based_community_detection(gr)
#'   e_gr<-add_system_group(e_gr, zendesk_group_filename)
#'   sn_in_cytpscape(e_gr)
#' }
