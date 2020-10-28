#' @title Get parcel-based disconnection
#' @description This function computes parcel-based direct disconnection measures using an
#' MNI-registered lesion file and an MNI-registered brain parcellation.
#' @param cfg a pre-made cfg structure (as list object).
#'
#' @importFrom R.matlab readMat
#'
#' @return A .trk.gz file. This contains all of the streamlines that intersected the lesion and can be viewed using e.g. DSI_Studio;
#' an .RData file with the suffix .connectivity.RData. This contains both the structural disconnection matrix (connectivity) and parcel names (name);
#' an .RData file with the suffix .network_measures.RData, which contains various graph measures for the SC matrix;
#' a .txt file with the suffix .connectogram.txt. This file contains a connectomgram that can be viewed on http://mkweb.bcgsc.ca/tableviewer/visualize/ by checking the two size options in step 2A (col with row size, row with col size);
#' a .mat file with the suffix _percent_parcel_SDC.mat. This file contains a disconnection adjacency matrix where each cell quantifies the % disconnection for each edge in the SC atlas;
#' a .mat file with the suffix _percent_parcel_spared_SC.mat. This file contains a spared connection adjacency matrix where each cell quantifies the % of each connection spared by the lesion for each edge in the SC atlas (i.e. it is equal to the difference between the Atlas SC matrix and the percent SDC matrix);
#' a .node file with the suffix _percent_parcel_SDC.node. This file contains the node information for external connectome viewers (e.g. MRIcroGL). Node sizes are proportional to the number of affected connections. Node colors can be pre-assigned in the .cfg file (cfg.node_color), but if not, they will be proportional to the amount of disconnection sustained analogous to node size;
#' a .edge file with the suffix _percent_parcel_SDC.edge. This file contains the percent SDC matrix in a format that can be loaded into external viewers (e.g. MRICroGL);
#' a .node file with the suffix _percent_parcel_spared_SC.node. This is analogous to (7), but for the spared SC matrix;
#' a .edge file with the suffix _percent_parcel_spared_SC.edge. This is analogous to (8), but for the spared SC matrix;
#' a .tdi.nii.gz file named the same way as the .trk.gz file. This contains a nifti image volume with track density imaging (TDI) values from the .trk.gz file at each voxel. It is essentially a way of converting the .trk.gz file into voxel space. Higher values indicate higher streamline densities at each grid element (voxel);
#' a .nii file with the suffix _percent_tdi.nii. For each voxel, values correspond the % reduction in streamline density relative to the atlas when accounting for the effects of the lesion.

#'
#' @examples \dontrun{
#'
#' }
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

  out_file = paste0(pd.path,cfg$pat_id,"_",cfg$file_suffix,'.trk.gz')

  system(paste0('! ',cfg$dsi_path,' --action=ana --source=',cfg$source_path,'/HCP842_1mm.fib.gz",
                " --tract=',cfg$source_path,'/all_tracts_1mm.trk.gz'," --roi=",cfg$lesion_path,
                ' --output=',out_file,' --connectivity=',cfg$parcel_path,' --connectivity_type=',
                cfg$con_type,' --connectivity_threshold=0 --export=tdi'))

  matfile=paste0(pd.path,"/",list.files(pd.path,pattern="connectivity\\.mat$"))
  mat=readMat(matfile)
  connectivity=mat$connectivity
  name=strsplit(intToUtf8(mat$name),"\n")[[1]]
  atlas=mat$atlas
  file.remove(matfile); rm(mat)

  netfile=paste0(pd.path,"/",list.files(pd.path,pattern="network_measures\\.txt$"))
  global=read.table(netfile,sep="\t")[1:27,]
  colnames(global)=c("Measure","Value")
  local=read.table(netfile,sep="\t",skip=27,header=T)[,-137]
  colnames(local)[1]="Measure"

  save(global,local,file=gsub("\\.txt","\\.RData",netfile))
  save(connectivity,name,atlas,file=gsub("\\.mat","\\.RData",matfile))

  pat_con=connectivity; rm(connectivity)




  cat('Finished creating atlas files')
}
