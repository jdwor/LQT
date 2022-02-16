#' @title Get parcel-based structural connectivity
#' @description This function computes parcel-based direct structural connectivity measures using an MNI-registered
#' brain parcellation and the curated HCP-842 structural connectome template described in Yeh et al., (2018 - NeuroImage)
#' @param cfg a pre-made cfg structure (as list object).
#'
#' @importFrom R.matlab readMat
#' @importFrom utils read.table
#'
#' @return An .RData file with the suffix .connectivity.RData, which contains the structural connection matrix (connectivity);
#' an .RData file with the suffix .network_measures.RData, which contains various graph measures for the SC matrix;
#' a .txt file with the suffix .connectogram.txt, which contains a connectogram that can be viewed on http://mkweb.bcgsc.ca/tableviewer/visualize/ by checking the two size options in step 2A (col with row size, row with col size);
#' a .tdi.nii.gz file named the same way as the .trk.gz file, which contains a nifti image with track density imaging (TDI) values at each voxel. It is essentially a way of converting the .trk.gz file into voxel space. Higher values indicate higher streamline densities at each voxel.
#'
#' @export

get_parcel_atlas<-function(cfg){
  if(is.null(cfg$pat_id)){
    at.path=paste0(cfg[[1]]$out_path,"/Atlas")
  }else{
    at.path=paste0(cfg$out_path,"/Atlas")
  }
  if(!dir.exists(at.path)){
    dir.create(at.path)
  }

  out_file = at.path

  out=tryCatch({
    suppressWarnings(system(paste0('! ',cfg$dsi_path,' --action=ana --source=',cfg$source_path,'/HCP842_1mm.fib.gz',
                                   ' --tract=',cfg$source_path,'/all_tracts_1mm.trk.gz',' --output=',out_file,
                                   ' --connectivity=',cfg$parcel_path,' --connectivity_type=',cfg$con_type,
                                   ' --connectivity_threshold=0 --export=tdi'),intern=T))
  },error=function(e){
    suppressWarnings(shell(paste0(cfg$dsi_path,' --action=ana --source=',cfg$source_path,'/HCP842_1mm.fib.gz',
                                   ' --tract=',cfg$source_path,'/all_tracts_1mm.trk.gz',' --output=',out_file,
                                   ' --connectivity=',cfg$parcel_path,' --connectivity_type=',cfg$con_type,
                                   ' --connectivity_threshold=0 --export=tdi'),intern=T))
  })

  oldnames=list.files(at.path, pattern="all_tracts_1mm",full.names=T)
  newnames=gsub("all_tracts_1mm",paste0('atlas_',cfg$file_suffix),oldnames)
  file.rename(oldnames,newnames)

  matfile=paste0(at.path,"/",list.files(at.path,pattern="connectivity\\.mat$"))
  mat=readMat(matfile)
  connectivity=mat$connectivity
  rownames(connectivity)=cfg$node_label
  colnames(connectivity)=cfg$node_label
  node_label=cfg$node_label
  node_group=cfg$node_group
  file.remove(matfile); rm(mat)

  netfile=paste0(at.path,"/",list.files(at.path,pattern="network_measures\\.txt$"))
  global=read.table(netfile,sep="\t")[1:27,]
  colnames(global)=c("Measure","Value")
  local=read.table(netfile,sep="\t",skip=27,header=T)[,-(length(node_label)+2)]
  colnames(local)[1]="Measure"
  colnames(local)[-1]=cfg$node_label
  file.remove(netfile)

  save(global,local,file=paste0(at.path,"/atlas_",cfg$file_suffix,"_network_measures.RData"))
  save(connectivity,node_label,node_group,
       file=paste0(at.path,"/atlas_",cfg$file_suffix,"_connectivity.RData"))

  cat('Finished creating atlas files.\n')
}
