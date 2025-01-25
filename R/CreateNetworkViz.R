#' Create Network Visualizations
#' 
#' This function creates a network visualization based on pathway data and machine learning results.
#' It generates a network of elements (genes, proteins, etc.) with edges representing interactions or relationships.
#' 
#' @param dat Data frame containing the data to visualize.
#' @param MLN_results Results from the multilevel network recovery model.
#' @param num_sets Number of sets (pathways) to visualize.
#' @param set_names A character vector containing names of the sets (pathways) to be used in the visualization.
#' 
#' @return A ggplot object with the network visualization.
#' @import class ggplot2 GGally network sna grpreg ggnetwork
#'
#' @export

CreateNetworkViz = function(dat, MLN_results, num_sets, set_names = c("Alanine and aspartate metabolism", "ATP synthesis","c22_U133_probes(user defined)","c23_U133_probes(user defined)","JAK-Stat Signaling Pathway","Oxidative phosphorylation","OXPHOS_HG-U133A_probes(user defined)","Parkinson's disease","Ubiquinone biosynthesis","gamma-Hexachlorocyclohexane degradation")){
  
  group_list = paste(c(1:num_sets))
  pwy_dfs = pathway_creator(dat[2:dim(dat)[2]], group_list, num_sets)
  og_df = dat[-c((length(y))),]
  
  total_genes = Reduce("+",lapply(pwy_dfs, function(x) {x <- dim(x)[2]}))
  rm(c)
  node_mat = data.frame(matrix(NA, nrow = total_genes, ncol = 3))
  node_mat[1:total_genes,1] = do.call(c, lapply(pwy_dfs, function(x) {x <- colnames(x)}))
  colnames(node_mat) = c("element", "set", "status")
  rownames(node_mat) <- as.character(node_mat$element)
  cntr = 0
  for(j in 1:num_pw){
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
  
  # pathway affiliation
  x = data.frame(element = network.vertex.names(net))
  x = merge(x, node_mat, by = "element", sort = FALSE)$set
  net %v% "set" = as.character(x)
  
  
  viz_net_obj = ggnetwork(net, arrow.gap = 0.02, layout = "kamadakawai")
  viz_net_obj_mod = viz_net_obj
  
  viz_net_obj_mod$x = viz_net_obj_mod$x+num_pw*(cos((as.numeric(viz_net_obj_mod$pathway)*(2*pi/num_pw))))/2.25
  viz_net_obj_mod$y = viz_net_obj_mod$y+num_pw*(sin((as.numeric(viz_net_obj_mod$pathway)*(2*pi/num_pw))))/2.25
  viz_net_obj_mod$xend = viz_net_obj_mod$xend+num_pw*(cos((as.numeric(viz_net_obj_mod$pathway)*(2*pi/num_pw))))/2.25
  viz_net_obj_mod$yend = viz_net_obj_mod$yend+num_pw*(sin((as.numeric(viz_net_obj_mod$pathway)*(2*pi/num_pw))))/2.25
  
  
  # looking up names by target coordinates
  viz_net_obj_mod$target.names = knn(train = viz_net_obj_mod[,c('x','y')], test = viz_net_obj_mod[,c('xend','yend')], cl = viz_net_obj_mod[,c('vertex.names')], k=1)
  
  elim_list = list()
  
  viz_net_obj_mod_rest = viz_net_obj_mod[!row(viz_net_obj_mod)[,1] %in% elim_list,]
  
  centroid_df = viz_net_obj_mod_rest %>%
    group_by(pathway) %>%
    dplyr::summarise(
      count = n(),
      x = (max(x,na.rm=TRUE)-min(x,na.rm=TRUE))/2+min(x,na.rm=TRUE),
      xend = (max(x,na.rm=TRUE)-min(x,na.rm=TRUE))/2+min(x,na.rm=TRUE),
      y = (max(y,na.rm=TRUE)-min(y,na.rm=TRUE))/2+min(y,na.rm=TRUE),
      yend = (max(y,na.rm=TRUE)-min(y,na.rm=TRUE))/2+min(y,na.rm=TRUE)
    )
  
  centroid_df = centroid_df[order(as.numeric(centroid_df$pathway)),]
  centroid_df$paired = gam_mod
  
  centroid_df_lines = centroid_df[centroid_df$paired==1,]
  idx <- combinations(sum(gam_mod), 2)
  idx = matrix(rbind(idx[,1], idx[,2]),ncol=1)
  centroid_df_lines <- centroid_df_lines[idx,]
  centroid_df_lines$paired = rep(1:choose(sum(gam_mod), 2), each = 2)
  
  group_list_names = set_names
  order_vec = as.numeric(viz_net_obj_mod_rest$pathway)
  group_list_names = rep(group_list_names, times = as.numeric(table(order_vec)))
  viz_net_obj_mod_rest2 = viz_net_obj_mod_rest
  viz_net_obj_mod_rest2$pathway = group_list_names
  
  ggplot(viz_net_obj_mod_rest2, 
         aes(x, y, xend = xend,yend = yend)) +
    geom_point() +
    geom_point(data = centroid_df,  shape = 21, col= "black", size=((5)+0.7)^2, stroke = 2)+
    geom_line(data = centroid_df_lines, aes(group = paired)) +
    geom_edges(aes(color = pathway), alpha = 0.1) +
    geom_nodes(aes(color = pathway), size = 0.75) +
    
    scale_color_brewer("Pathway", palette = "Spectral")+
    xlim(-num_pw*0.6, num_pw*0.6)+
    ylim(-num_pw*0.6, num_pw*0.6)+
    xlab("")+
    ylab("")+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
}
