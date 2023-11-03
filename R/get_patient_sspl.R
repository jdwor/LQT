#' @title Get structural shortest path length (SSPL)
#' @description This function computes increases in the structural shortest path length (SSPL) between brain regions
#' associated with focal lesions. See Griffis et al. (2019, NeuroImage) for more details.
#' @param cfg a pre-made cfg structure (as list object).
#' @param cores an integer value that indicates how many parallel cores the function should be run on.
#'
#' @importFrom igraph distances graph_from_adjacency_matrix V
#' @importFrom pbmcapply pbmclapply
#'
#' @return A file with the suffix _SSPL_matrices.RData -- this contains the raw SSPL matrix (sspl_matrix), the lesion-based increases in SSPLs relative to the atlas (both direct and indirect disconnections; delta_ssplL_matrix), and the indirect-only delta matrix (idc_matrix);
#' a file with the suffix _delta_SSPL.node -- this contains the node information for loading the delta SSPL matrix into an external viewer. Node sizes are proportional to the total SSPL increases for each node (relative to the max increase);
#' a file with the suffix _delta_SSPL.edge -- this contains the edge information for loading the delta SSPL matrix into an external viewer;
#' a file with the suffix _indirect_SDC.node -- same as (3), but for the indirect SDC matrix;
#' a file with the suffix _indirect_SDC.edge -- same as (4), but for the indirect SDC matrix.
#'
#' Note that if the lesion leads to the SSPL between a region pair to be undefined (i.e. = Inf), then this value is replaced
#' by a value equal to 1 plus the maximum non-Inf value in the atlas SSPL matrix.
#'
#' @export

get_patient_sspl<-function(cfg, cores=1, verbose=T){
  if(is.null(cfg$pat_id)){

    at.path=paste0(cfg[[1]]$out_path,"/Atlas")
    if(!dir.exists(at.path)){
      stop('Atlas folder does not exist in output directory; cannot compute changes in SSPLs until atlas SSPL matrix has been created. Please create atlas SSPL matrix and retry.');
    }
    cat('Computing patient SSPL matrices\n')
    out = pbmclapply(cfg, get_patient_sspl, verbose=F, mc.cores=cores)
    cat('Finished creating patient SSPL and delta-SSPL matrices\n')

  }else{
    at.path=paste0(cfg$out_path,"/Atlas")
    if(!dir.exists(at.path)){
      stop('Atlas folder does not exist in output directory; cannot compute changes in SSPLs until atlas SSPL matrix has been created. Please create atlas SSPL matrix and retry.');
    }

    # load atlas SSPL matrix
    load(paste0(at.path,'/atlas_',cfg$file_suffix,'_SSPL_matrix.RData'))
    atlas_sspl = sspl_matrix; rm(sspl_matrix)

    # Load patient SC and SDC matrices
    pat.path=paste0(cfg$out_path,"/",cfg$pat_id)
    pd.path=paste0(pat.path,"/Parcel_Disconnection")
    if(!dir.exists(pat.path)){
      stop('Patient folder does not exist in output directory; cannot compute changes in SSPLs until patient disconnectivity has been calculated.');
    }
    if(!dir.exists(pd.path)){
      stop('Patient Parcel Disconnection folder does not exist; cannot compute changes in SSPLs until patient disconnectivity has been calculated.');
    }

    load(paste0(pd.path,"/",cfg$pat_id,'_',cfg$file_suffix,'_percent_parcel_mats.RData'))

    ps.path=paste0(pat.path,"/Parcel_SSPL")
    if(!dir.exists(ps.path)){
      dir.create(ps.path)
    }

    # get SSPLs for patient
    if(is.null(cfg$sspl_spared_thresh)){
      cfg$sspl_spared_thresh=10
    }else if(cfg$sspl_spared_thresh>100 | cfg$sspl_spared_thresh<1){
      stop('Error: SSPL threshold out of bounds. Please check and try again')
    }

    graph = graph_from_adjacency_matrix(pct_spared_sc_matrix>=cfg$sspl_spared_thresh,mode="undirected")
    sspl_matrix = distances(graph,v=V(graph),to=V(graph)) # SSPL for end atlas
    diag(sspl_matrix)=NA # set diagonal to NA
    sspl_matrix[is.na(atlas_sspl)] = NA

    # get delta SSPLs for patient
    sspl_matrix[is.infinite(sspl_matrix)] = max(as.vector(atlas_sspl),na.rm=T)+1
    delta_sspl_matrix = sspl_matrix-atlas_sspl
    delta_sspl_matrix[is.na(delta_sspl_matrix)]=0

    # get parcel coordinates if not supplied
    if(is.null(cfg$parcel_coords)){
      cat('Warning: Parcel coordinates not provided, attempting to estimate coordinates from parcellation image\n');
      cfg$parcel_coords = util_get_coords(cfg)
    }
    NodePos = cfg$parcel_coords

    # write out .edge file
    write(round(t(delta_sspl_matrix),4),
          paste0(ps.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                 "_delta_SSPL.edge"),
          ncolumns=ncol(delta_sspl_matrix),sep="\t")

    # write out .node files
    NodeSize = apply(sspl_matrix,2,sum,na.rm=T)/apply(atlas_sspl,2,sum,na.rm=T) # size nodes according to % maximum # of affected connections
    NodeSize[is.na(NodeSize)]=0
    if(is.null(cfg$node_color)){
      NodeColor = NodeSize
    }else if(length(cfg$node_color)!=length(NodeSize)){
      NodeColor = NodeSize
    }else{
      NodeColor = cfg$node_color
    }
    if(is.null(cfg$node_label)){
      NodeLabel = as.character(1:nrow(delta_sspl_matrix))
    }else if(length(cfg$node_label)!=length(NodeSize)){
      NodeLabel = as.character(1:nrow(delta_sspl_matrix))
    }else{
      NodeLabel = cfg$node_label
    }

    # construct output filename and save
    OutputFile = paste0(ps.path,"/",cfg$pat_id,"_",cfg$file_suffix,"_delta_SSPL.node")
    C = cbind(NodePos, NodeColor, NodeSize, NodeLabel)
    write(t(C), OutputFile, sep="\t")

    if(verbose==T){cat('Finished creating patient SSPL and delta-SSPL matrices\n')}

    ### output direct and indirect structural disconnection matrices
    idc_matrix_mask = (delta_sspl_matrix > 0) - (pct_sdc_matrix > 0)
    idc_matrix = delta_sspl_matrix
    idc_matrix[idc_matrix_mask==0]=0
    idc_matrix[is.na(idc_matrix)]=0

    rownames(sspl_matrix)=cfg$node_label
    colnames(sspl_matrix)=cfg$node_label
    rownames(delta_sspl_matrix)=cfg$node_label
    colnames(delta_sspl_matrix)=cfg$node_label
    rownames(idc_matrix)=cfg$node_label
    colnames(idc_matrix)=cfg$node_label
    node_label=cfg$node_label
    node_group=cfg$node_group

    save(sspl_matrix,delta_sspl_matrix,idc_matrix,node_label,node_group,
         file=paste0(ps.path,"/",cfg$pat_id,"_",cfg$file_suffix,"_SSPL_matrices.RData"))

    ### save parcel-level sspl averages
    
    parc.sspl=data.frame(Parcel=cfg$node_label,ParcelGroup=cfg$node_group)
    parc.sspl$Parcel=gsub(" \\(Right\\)","_R",parc.sspl$Parcel)
    parc.sspl$Parcel=gsub(" \\(Left\\)","_L",parc.sspl$Parcel)
    parc.sspl$Parcel=gsub("Lenticular nucleus, p",
                         "Lent_Nuc_P",parc.sspl$Parcel)
    parc.sspl$Parcel=gsub("RH_","R_",parc.sspl$Parcel)
    parc.sspl$Parcel=gsub("LH_","L_",parc.sspl$Parcel)

    ps.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_SSPL")
    load(paste0(ps.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                "_SSPL_matrices.RData"))
    perc_d_sspl=delta_sspl_matrix/(sspl_matrix-delta_sspl_matrix)
    perc_d_sspl[is.na(perc_d_sspl)]=0

    parc.sspl$ParcelAvgSSPL=apply(perc_d_sspl,1,mean)
    write.csv(parc.sspl,paste0(ps.path,"/",cfg$pat_id,"_",
                                 cfg$file_suffix,"_parcel_average_sspl.csv"))

    # write out .edge file
    write(round(t(idc_matrix),4),
          paste0(ps.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                 "_indirect_SDC_edge"),
          ncolumns=ncol(idc_matrix),sep="\t")

    # construct output filename and save
    OutputFile = paste0(ps.path,"/",cfg$pat_id,"_",cfg$file_suffix,"_indirect_SDC.node")
    NodeSize = apply(idc_matrix,2,sum,na.rm=T)/(max(apply(idc_matrix,2,sum,na.rm=T))+.01) # size nodes according to % maximum # of affected connections in matrix
    if(is.null(cfg$node_color)){
      NodeColor = NodeSize
    }else if(length(cfg$node_color)!=length(NodeSize)){
      NodeColor = NodeSize
    }else{
      NodeColor = cfg$node_color
    }
    C = cbind(NodePos, NodeColor, NodeSize, NodeLabel)
    write(t(C), OutputFile, sep="\t")

  }

}


