#' @title Get parcel-based structural connectivity
#' @description This function computes parcel-based direct structural connectivity measures using an MNI-registered
#' brain parcellation and the curated HCP-842 structural connectome template described in Yeh et al., (2018 - NeuroImage)
#' @param cfg a pre-made cfg structure (as list object).
#' @param saveout is a logical value that indicates whether the user would like to save output files to cfg$out_path
#'
#' @importFrom neurobase readnii writenii
#' @importFrom R.matlab readMat
#'
#' @return A data.frame giving the percent of streamlines in each tract that were disconnected by the lesions (column "Discon")
#' and the associated tract names (column "Tract").
#' @examples \dontrun{
#'
#' }
#' @export

get_parcel_atlas<-function(cfg,saveout=TRUE){
  at.path=paste0(cfg$source_path,"/Atlas")
  if(!dir.exists(at.path)){
    dir.create(at.path)
  }

  out_file = paste0(at.path,'atlas_',cfg$file_suffix,'.trk.gz')

  system(paste0('! ',cfg$dsi_path,' --action=ana --source=',cfg$source_path,'/HCP842_1mm.fib.gz",
                " --tract=',cfg$source_path,'/all_tracts_1mm.trk.gz',' --output=',out_file,' --connectivity=',cfg$parcel_path,
                ' --connectivity_type=',cfg$con_type,' --connectivity_threshold=0 --export=tdi'))

  cat('Finished creating atlas files')
}
