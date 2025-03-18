#' Create Network Visualizations
#'
#' This function creates a network visualization based on set data and machine learning results.
#' It generates a network of elements (genes, proteins, etc.) with edges representing interactions or relationships.
#'
#' @param dat Data frame containing the data to visualize.
#' @param gamma set Results from the multilevel network recovery model.
#' @param xi element Results from the multilevel network recovery model.
#' @param num_sets Number of sets (pathways) to visualize.
#' @param set_names A character vector containing names of the sets (pathways) to be used in the visualization.
#'
#' @return A ggplot object with the network visualization.
#' @import GGally
#' @import GIGrvg
#' @import class
#' @import dplyr
#' @import expm
#' @import fMultivar
#' @import ggnetwork
#' @import ggplot2
#' @import gtools
#' @import infotheo
#' @import invgamma
#' @import laGP
#' @import mvtnorm
#' @import network
#' @import rpart
#' @import sna
#' @import tidyr
#' @import tidyverse
#' @importFrom Matrix forceSymmetric
#' @importFrom ClusterR GMM
#' @importFrom ClusterR predict_GMM
#' @importFrom ald rALD dALD
#' @importFrom plgp covar
#' @importFrom plgp covar.sep
#' @importFrom utils combn
#' @export
CreateNetworkViz = function(dat, gamma, xi, num_sets, set_names = NULL){

  y = dat[1:(dim(dat)[1]-1),1]
  gam_mod = gamma
  MLN_results = xi

  group_list = paste(c(1:num_sets))
  pwy_dfs = pathway_creator(dat, num_sets)
  og_df = dat[-c((length(y))),]

  total_genes = Reduce("+",lapply(pwy_dfs, function(x) {x <- dim(x)[2]}))
  node_mat = data.frame(matrix(NA, nrow = total_genes, ncol = 3))
  node_mat[1:total_genes,1] = do.call(c, lapply(pwy_dfs, function(x) {x <- colnames(x)}))
  colnames(node_mat) = c("element", "set", "status")
  rownames(node_mat) <- as.character(node_mat$element)
  cntr = 0
  for(j in 1:num_sets){
    node_mat[c(colnames(pwy_dfs[[j]])),2] = j
    node_mat[c(colnames(pwy_dfs[[j]])),3] = as.numeric(MLN_results[cntr+1:dim(pwy_dfs[[j]])[2]])
    cntr = cntr + dim(pwy_dfs[[j]])[2]
  }

  edge_list = singular_pts = data.frame(V1 = character(), V2 = character())
  for(np in 1:num_sets){
    df_plc = node_mat[((node_mat$set==np)&(node_mat$status==1)),]
    if(nrow(df_plc>1)){
      df_ST = t(combn(df_plc$element,2))
      edge_list = rbind(edge_list,as.data.frame(df_ST))
    } else if(nrow(df_plc)==1){
      df_SP = data.frame(V1 = df_plc$element, V2 = df_plc$element)
      singular_pts = rbind(singular_pts,df_SP)
    }
    df_plc = node_mat[((node_mat$set==np)&(node_mat$status==0)),]
    df_SP = data.frame(V1 = df_plc$element, V2 = df_plc$element)
    singular_pts = rbind(singular_pts,df_SP)
  }

  full_edge = rbind(edge_list,singular_pts)
  edge_mat = full_edge
  colnames(edge_mat) = c("Source", "Target")

  # network object
  net = network(edge_mat, directed = FALSE, loops = TRUE)

  # set affiliation
  x = data.frame(element = network.vertex.names(net))
  x = merge(x, node_mat, by = "element", sort = FALSE)$set
  net %v% "set" = as.character(x)


  viz_net_obj = ggnetwork(net, arrow.gap = 0.02, layout = "kamadakawai")
  viz_net_obj_mod = viz_net_obj

  viz_net_obj_mod$x = viz_net_obj_mod$x+num_sets*(cos((as.numeric(viz_net_obj_mod$set)*(2*pi/num_sets))))/2.25
  viz_net_obj_mod$y = viz_net_obj_mod$y+num_sets*(sin((as.numeric(viz_net_obj_mod$set)*(2*pi/num_sets))))/2.25
  viz_net_obj_mod$xend = viz_net_obj_mod$xend+num_sets*(cos((as.numeric(viz_net_obj_mod$set)*(2*pi/num_sets))))/2.25
  viz_net_obj_mod$yend = viz_net_obj_mod$yend+num_sets*(sin((as.numeric(viz_net_obj_mod$set)*(2*pi/num_sets))))/2.25

  # looking up names by target coordinates
  viz_net_obj_mod$target.names = knn(train = viz_net_obj_mod[,c('x','y')], test = viz_net_obj_mod[,c('xend','yend')], cl = viz_net_obj_mod[,c('vertex.names')], k=1)

  elim_list = list()

  viz_net_obj_mod_rest = viz_net_obj_mod[!row(viz_net_obj_mod)[,1] %in% elim_list,]

  radius_df = viz_net_obj_mod_rest %>%
    group_by(set) %>%
    mutate(
      x_center = mean(x, na.rm = TRUE),
      y_center = mean(y, na.rm = TRUE),
      distance = sqrt((x - x_center)^2 + (y - y_center)^2)
    ) %>%
    summarise(
      count = n(),
      x = mean(x, na.rm = TRUE),
      xend = mean(x, na.rm = TRUE),
      y = mean(y, na.rm = TRUE),
      yend = mean(y, na.rm = TRUE),
      radius = max(distance, na.rm = TRUE)  # Max distance from centroid
    )

  centroid_df = viz_net_obj_mod_rest %>%
    group_by(set) %>%
    dplyr::summarise(
      count = n(),
      x = (max(x,na.rm=TRUE)-min(x,na.rm=TRUE))/2+min(x,na.rm=TRUE),
      xend = (max(x,na.rm=TRUE)-min(x,na.rm=TRUE))/2+min(x,na.rm=TRUE),
      y = (max(y,na.rm=TRUE)-min(y,na.rm=TRUE))/2+min(y,na.rm=TRUE),
      yend = (max(y,na.rm=TRUE)-min(y,na.rm=TRUE))/2+min(y,na.rm=TRUE),
      radius = max(sqrt((x - mean(x, na.rm=TRUE))^2 + (y - mean(y, na.rm=TRUE))^2)) # Max distance from center
    )

  centroid_df$radius = radius_df$radius
  centroid_df = centroid_df[order(as.numeric(centroid_df$set)),]
  centroid_df$paired = gam_mod

  if(sum(gam_mod)>1){
    centroid_df_lines = centroid_df[centroid_df$paired==1,]
    idx <- combinations(sum(gam_mod), 2)
    idx = matrix(rbind(idx[,1], idx[,2]),ncol=1)
    centroid_df_lines <- centroid_df_lines[idx,]
    centroid_df_lines$paired = rep(1:choose(sum(gam_mod), 2), each = 2)

    group_list_names = set_names
    order_vec = as.numeric(viz_net_obj_mod_rest$set)
    group_list_names = rep(group_list_names, times = as.numeric(table(order_vec)))
    viz_net_obj_mod_rest2 = viz_net_obj_mod_rest
    viz_net_obj_mod_rest2$set = group_list_names

    ggplot(viz_net_obj_mod_rest2,
           aes(x, y, xend = xend,yend = yend)) +
      geom_point() +
      geom_point(data = centroid_df,  shape = 21, col= "black", size=centroid_df$radius * 75, stroke = 2)+
      geom_line(data = centroid_df_lines, aes(group = paired)) +
      geom_edges(aes(color = set), alpha = 0.1) +
      geom_nodes(aes(color = set), size = 0.75) +

      scale_color_brewer("set", palette = "Spectral")+
      xlim(-num_sets*0.6, num_sets*0.6)+
      ylim(-num_sets*0.6, num_sets*0.6)+
      xlab("")+
      ylab("")+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank())

  } else {

    group_list_names = set_names
    order_vec = as.numeric(viz_net_obj_mod_rest$set)
    group_list_names = rep(group_list_names, times = as.numeric(table(order_vec)))
    # viz_net_obj_mod_rest2 = viz_net_obj_mod_rest
    # viz_net_obj_mod_rest2$set = group_list_names
    viz_net_obj_mod_rest2 <- merge(
      viz_net_obj_mod_rest,
      node_mat[, c("element", "set")],
      by.x = "vertex.names",
      by.y = "element",
      all.x = TRUE
    )

    # If `set.y` is the correct one (from node_mat), drop `set.x`
    viz_net_obj_mod_rest2$set <- viz_net_obj_mod_rest2$set.y
    viz_net_obj_mod_rest2$set.x <- NULL
    viz_net_obj_mod_rest2$set.y <- NULL

    ggplot(viz_net_obj_mod_rest2,
           aes(x, y, xend = xend,yend = yend)) +
      geom_point() +
      geom_point(data = centroid_df,  shape = 21, col= "black", size=centroid_df$radius * 75, stroke = 2)+
      geom_nodes(aes(color = as.factor(set)), size = 0.75) +
      geom_edges(aes(color = as.factor(set)), alpha = 0.1) +

      scale_color_brewer("set", palette = "Spectral")+
      xlim(-num_sets*0.6, num_sets*0.6)+
      ylim(-num_sets*0.6, num_sets*0.6)+
      xlab("")+
      ylab("")+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
  }

}
