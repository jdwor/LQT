#' @title Compute the structural shortest path length (SSPL) for the SC atlas
#' @description This function uses the atlas SC matrix to compute an atlas SSPL matrix
#' containing SSPLs between each pair of brain regions. It assumes that you have already
#' run get_parcel_atlas to obtain the atlas SC matrix.
#' @param cfg a pre-made cfg structure (as list object).
#'
#' @importFrom igraph distances graph_from_adjacency_matrix V
#'
#' @return An .RData file with the suffix _sspl_matrix.RData, which contains the raw SSPL
#' matrix computed from the binarized atlas SC matrix.
#'
#' @export

get_atlas_sspl<-function(cfg){
  at.path=paste0(cfg$out_path,"/Atlas")
  if(!dir.exists(at.path)){
    stop('Atlas folder does not exist in output directory; cannot compute atlas SSPLs until atlas SC matrix has been created. Please create atlas SC matrix and retry.');
  }

  # find and load SC matrix file
  load(paste0(at.path,'/atlas_',cfg$file_suffix,'_connectivity.RData'))

  # get SSPLs for atlas
  graph = graph_from_adjacency_matrix(connectivity,mode="undirected")
  sspl_matrix = distances(graph,v=V(graph),to=V(graph)) # SSPL for end atlas
  diag(sspl_matrix)=NA # set diagonal to NA
  sspl_matrix[is.infinite(sspl_matrix)] = NA # set any Inf values to NA
  rownames(sspl_matrix)=cfg$node_label
  colnames(sspl_matrix)=cfg$node_label
  node_group=cfg$node_group

  out_file = paste0(at.path,'/atlas_',cfg$file_suffix,'_SSPL_matrix.RData')
  save(sspl_matrix,node_group,file=out_file)

  cat('Finished calculating atlas SSPL.\n')
}
