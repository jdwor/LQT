#' @title Get parcel damage
#' @description % This function uses an MNI-registered lesion file and an MNI-registered brain parcellation
#' to estimate the amount of damage sustained by each brain region.
#' @param cfg a pre-made cfg structure (as list object).
#' @param saveout is a logical value that indicates whether the user would like to save output files to cfg$out_path
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel.
#' @param cores if parallel = TRUE, cores is an integer value that indicates how many cores
#' the function should be run on.
#' @param pb is a logical value that indicates whether or not a progress bar will be shown during analysis.
#'
#' @importFrom neurobase readnii writenii
#'
#' @return A list containing
#' 1) damage.map (a nifti image in which voxel values correspond to the % of the associated parcel damaged by lesion),
#' and 2) damage.vec (a vector with % damage values for each parcel in the brain parcellation).
#' @examples \dontrun{
#' mmdt.obj = get.mmdt.obj(masks = masks, modal1 = t1s, modal2 = flairs,
#'                         ids = ids, groups = groups)
#'
#' mmdt.results = mmdt(mmdt.obj, mc.adjust = c("BH", "maxt"), nperm = 500,
#'                     parallel = TRUE, cores = 4, pb = TRUE)}
#' @export

get_parcel_damage<-function(cfg,saveout=TRUE,parallel=TRUE,cores=2,pb=TRUE){

  lesions = readnii(cfg$lesion_path)
  parcels = readnii(cfg$parcel_path)

  # Check if lesion and parcels are the same dimensions
  if(!identical(dim(lesion.img),dim(parcels.img))){
    stop('Lesion and parcel volumes have different dimensions. Please make sure that both volumes are in the same space and have identical dimensions.')
  }
  else{
    cat('Lesion and parcel volumes have identical dimensions; computing parcel damage measures.')
  }

  # Get parcel damage measures
  pcd_vol = parcels;
  pcd_vol[pcd_vol>0]=0
  num_parcels = length(unique(parcels[parcels>0]))
  pcd_vect = rep(0, num_parcels)
  for(i in 1:num_parcels){
    # loop through parcels

    cat('Progress:',i,'of',num_parcels,'parcels evaluated\n')
    pcd_vect[i] = mean(lesions[parcels==i])*100 # get percent of parcel damaged
    pcd_vol[parcels==i] = pcd_vect[i] # put into volume
  }

  cat('Finished computing parcel damage measures')

  # save parcel damage measures

  if(saveout==T){
    pat.path=paste0(cfg$out_path,"/",cfg$pat_id)
    if(!dir.exists(cfg$out_path)){
      dir.create(cfg$out_path)
    }
    if(!dir.exists(pat.path)){
      dir.create(pat.path)
    }

    pd.path=paste0(pat.path,"/Parcel_Damage")
    if(!dir.exists(pd.path)){
      dir.create(pd.path)
    }

    writenii(lesion, paste0(pat.path,"/",cfg$pat_id,"_lesion.nii.gz"))
    writenii(pcd_vol, paste0(pat.path,"/",cfg$pat_id,"_",
                             cfg$file_suffix,"_percent_damage.nii.gz"))
    write.csv(pcd_vect,paste0(pat.path,"/",cfg$pat_id,"_",
                              cfg$file_suffix,"_percent_damage.csv"))
  }

  return(list(damage.map=pcd_vol, damage.vec=pcd_vect))
}
