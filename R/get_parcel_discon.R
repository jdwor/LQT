#' @title Get parcel-based disconnection
#' @description This function computes parcel-based direct disconnection measures using an
#' MNI-registered lesion file and an MNI-registered brain parcellation.
#' @param cfg a pre-made cfg structure (as list object).
#'
#' @importFrom R.matlab readMat
#' @importFrom neurobase readnii writenii
#' @importFrom fslr fslsmooth
#'
#' @return An .RData file with the suffix .connectivity.RData. This contains the structural disconnection matrix (connectivity);
#' an .RData file with the suffix .network_measures.RData, which contains various graph measures for the SC matrix;
#' an .RData file with the suffix _percent_parcel_mats.RData. This file contains a disconnection adjacency matrix (pct_sdc_matrix) and a spared connection adjacency matrix (pct_spared_sc_matrix);
#' a .txt file with the suffix .connectogram.txt. This file contains a connectogram that can be viewed on http://mkweb.bcgsc.ca/tableviewer/visualize/ by checking the two size options in step 2A (col with row size, row with col size);
#' a .node file with the suffix _percent_parcel_SDC.node. This file contains the node information for external connectome viewers (e.g. MRIcroGL). Node sizes are proportional to the number of affected connections. Node colors can be pre-assigned in the .cfg file (cfg.node_color), but if not, they will be proportional to the amount of disconnection sustained analogous to node size;
#' a .edge file with the suffix _percent_parcel_SDC.edge. This file contains the percent SDC matrix in a format that can be loaded into external viewers (e.g. MRICroGL);
#' a .node file with the suffix _percent_parcel_spared_SC.node. This is analogous to (7), but for the spared SC matrix;
#' a .edge file with the suffix _percent_parcel_spared_SC.edge. This is analogous to (8), but for the spared SC matrix;
#' a .tdi.nii.gz file named the same way as the .trk.gz file. This contains a nifti image volume with track density imaging (TDI) values from the .trk.gz file at each voxel. It is essentially a way of converting the .trk.gz file into voxel space. Higher values indicate higher streamline densities at each grid element (voxel);
#' a .nii file with the suffix _percent_tdi.nii.gz. For each voxel, values correspond the % reduction in streamline density relative to the atlas when accounting for the effects of the lesion.
#'
#' @export

get_parcel_discon<-function(cfg){
  at.path=paste0(cfg$out_path,"/Atlas")
  if(!dir.exists(at.path)){
    stop('Atlas folder does not exist in output directory; cannot convert streamline counts or TDI values to % disconnection. Please run util_get_parcel_atlas before attempting to create patient disconnection measures.');
  }
  pat.path=paste0(cfg$out_path,"/",cfg$pat_id)
  if(!dir.exists(pat.path)){
    dir.create(pat.path)
  }
  pd.path=paste0(pat.path,"/Parcel_Disconnection")
  dm.path=paste0(pat.path,"/Disconnection_Maps")
  if(!dir.exists(pd.path)){
    dir.create(pd.path)
  }
  if(!dir.exists(dm.path)){
    dir.create(dm.path)
  }

  out_file = paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,'.trk.gz')

  # create disconnection matrix, .trk fle, and TDI map
  out=suppressWarnings(system(paste0('! ',cfg$dsi_path,' --action=ana --source=',cfg$source_path,'/HCP842_1mm.fib.gz',
                ' --tract=',cfg$source_path,'/all_tracts_1mm.trk.gz',' --roi=',cfg$lesion_path,
                ' --output=',out_file,' --connectivity=',cfg$parcel_path,' --connectivity_type=',
                cfg$con_type,' --connectivity_threshold=0 --export=tdi'),intern=T))

  # resave connectivity file (raw streamline counts)
  matfile=paste0(pd.path,"/",list.files(pd.path,pattern="connectivity\\.mat$"))
  mat=readMat(matfile)
  connectivity=mat$connectivity
  rownames(connectivity)=cfg$node_label
  colnames(connectivity)=cfg$node_label
  node_group=cfg$node_group
  name=strsplit(intToUtf8(mat$name),"\n")[[1]]
  atlas=mat$atlas
  file.remove(matfile); rm(mat)

  # resave network measures file (raw streamline counts)
  netfile=paste0(pd.path,"/",list.files(pd.path,pattern="network_measures\\.txt$"))
  global=read.table(netfile,sep="\t")[1:27,]
  colnames(global)=c("Measure","Value")
  local=read.table(netfile,sep="\t",skip=27,header=T)[,-137]
  colnames(local)[1]="Measure"
  colnames(local)[-1]=cfg$node_label
  file.remove(netfile)

  save(global,local,file=gsub("\\.txt","\\.RData",netfile))
  save(connectivity,node_group,file=gsub("\\.mat","\\.RData",matfile))

  # load atlas SC matrix
  pat_con=connectivity; rm(connectivity)
  load(paste0(at.path,"/",list.files(at.path,pattern="connectivity\\.RData")))
  atlas_con=connectivity; rm(connectivity)

  # convert patient matrix to % disconnection and save
  pct_sdc_matrix = 100*(pat_con/atlas_con)
  pct_sdc_matrix[is.na(pct_sdc_matrix)]=0
  pct_spared_sc_matrix = 100*((atlas_con-pat_con)/atlas_con)
  pct_spared_sc_matrix[is.na(pct_spared_sc_matrix)]=0
  pct_spared_sc_matrix[is.infinite(pct_spared_sc_matrix)]=0

  rownames(pct_sdc_matrix)=cfg$node_label
  rownames(pct_sdc_matrix)=cfg$node_label
  colnames(pct_spared_sc_matrix)=cfg$node_label
  colnames(pct_spared_sc_matrix)=cfg$node_label

  fname=paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,"_percent_parcel_mats.RData")
  save(pct_sdc_matrix, pct_spared_sc_matrix, file=fname)

  ### output .node and .edge files for external viewers (e.g. MRIcroGL)

  # get parcel coordinates if not supplied
  if(is.null(cfg$parcel_coords)){
    cat('Warning: Parcel coordinates not provided, attempting to estimate coordinates from parcellation image\n');
    cfg$parcel_coords = util_get_coords(cfg)
  }
  NodePos = cfg$parcel_coords

  # write out .edge files
  write(round(t(pct_sdc_matrix),4),
        paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
               "_percent_parcel_SDC.edge"),
        ncolumns=ncol(pct_sdc_matrix),sep="\t")
  write(round(t(pct_spared_sc_matrix),4),
        paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
               "_percent_parcel_spared_SC.edge"),
        ncolumns=ncol(pct_spared_sc_matrix),sep="\t")

  # write out .node files
  NodeSize = apply(pct_sdc_matrix,2,sum,na.rm=T)/apply((atlas_con>0)*100,2,sum,na.rm=T) # size nodes according to % maximum # of affected connections
  NodeSize[is.na(NodeSize)]=0
  NodeSize2 = apply(pct_spared_sc_matrix,2,sum,na.rm=T)/apply((atlas_con>0)*100,2,sum,na.rm=T)
  NodeSize2[is.na(NodeSize2)]=0
  if(is.null(cfg$node_color)){
    NodeColor = NodeSize
  }else if(length(cfg$node_color)!=length(NodeSize)){
    NodeColor = NodeSize
  }else{
    NodeColor = cfg$node_color
  }
  if(is.null(cfg$node_label)){
    NodeLabel = as.character(1:nrow(pct_sdc_matrix))
  }else if(length(cfg$node_label)!=length(NodeSize)){
    NodeLabel = as.character(1:nrow(pct_sdc_matrix))
  }else{
    NodeLabel = cfg$node_label
  }

  # construct output filename and save
  OutputFile = paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,"_percent_parcel_SDC.node")
  C = cbind(NodePos, NodeColor, NodeSize, NodeLabel)
  write(t(C), OutputFile, sep="\t")

  OutputFile = paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,"_percent_parcel_spared_SC.node")
  C = cbind(NodePos, NodeColor, NodeSize2, NodeLabel)
  write(t(C), OutputFile, sep="\t")

  ### Save TDI map
  # convert TDI to % disconnection
  con_file = list.files(pd.path, pattern="\\.tdi\\.nii\\.gz$")
  pat_con = readnii(paste0(pd.path,"/",con_file))

  atlas_file = list.files(at.path, pattern="\\.tdi\\.nii\\.gz$")
  atlas_con = readnii(paste0(at.path,"/",atlas_file))

  writenii(pat_con, paste0(dm.path,"/",con_file))

  pat_con = 100*(pat_con/atlas_con)
  pat_con@datatype = 16
  support = substr(cfg$source_path,1,nchar(cfg$source_path)-18)
  qa = readnii(paste0(support,'HCP842_QA.nii.gz'))

  pat_con[is.na(pat_con)]=0
  pat_con[qa==0]=0
  pat_con[pat_con < 1] = 0

  if(!is.null(cfg$smooth)){
    if(cfg$smooth>0){
      pat_con = fslsmooth(pat_con,sigma=cfg$smooth,retimg=T,verbose=F)
      pat_con[pat_con < 0.01] = 0
    }
  }

  writenii(pat_con, paste0(dm.path,"/",cfg$pat_id,"_",cfg$file_suffix,"_percent_tdi.nii.gz"))

  trk_file = list.files(pd.path, pattern="\\.trk\\.gz$")
  if(length(trk_file)>0){file.remove(paste0(pd.path,"/",trk_file))}
  file.remove(paste0(pd.path,"/",con_file))

  cat('Finished computing patient disconnection measures.\n');
}
