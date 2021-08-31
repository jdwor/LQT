#' @title Get parcel damage
#' @description This function uses an MNI-registered lesion file and an MNI-registered brain parcellation
#' to estimate the amount of damage sustained by each brain region.
#' @param cfg a pre-made cfg structure (as list object).
#' @param cores an integer value that indicates how many parallel cores the function should be run on.
#'
#' @importFrom neurobase readnii writenii
#' @importFrom utils write.csv
#' @importFrom pbmcapply pbmclapply
#' @importFrom RNiftyReg buildAffine applyTransform
#'
#' @return A .nii.gz file with the suffix _percent_damage.nii.gz. This file contains a damage.map in which voxel values correspond to the % of the associated parcel damaged by lesion;
#' a .csv file with the suffix _percent_damage.csv. This file contains a data frame with % damage values for each parcel in the given parcellation.
#'
#' @export

get_parcel_damage<-function(cfg, cores=1, verbose=T){
  if(is.null(cfg$pat_id)){

    cat('Computing parcel damage measures.\n')
    out = pbmclapply(cfg, get_parcel_damage, verbose=F, mc.cores=cores)
    cat('Finished computing parcel damage measures')

  }else{

    lesions_or = readnii(cfg$lesion_path)
    trans = diag(rbind(lesions_or@srow_x,lesions_or@srow_y,
                       lesions_or@srow_z)[,1:3])
    if(!sum(trans)%in%c(0,3)){
      affine=buildAffine(scales=trans,source=lesions_or)
      lesions=applyTransform(affine,lesions_or)
    }else{
      lesions=lesions_or
    }
    parcels = readnii(cfg$parcel_path)

    # Check if lesion and parcels are the same dimensions
    if(!identical(dim(lesions),dim(parcels))){
      stop('Lesion and parcel volumes have different dimensions. Please make sure that both volumes are in the same space and have identical dimensions.')
    }else{
      if(verbose==T){cat('Computing parcel damage measures.\n')}
    }

    # Get parcel damage measures
    pcd_vol = parcels;
    pcd_vol[pcd_vol>0]=0
    num_parcels = length(unique(parcels[parcels>0]))
    pcd_vect = rep(0, num_parcels)
    for(i in 1:num_parcels){
      # loop through parcels

      if(verbose==T){cat('Progress:',i,'of',num_parcels,'parcels evaluated\r')}
      pcd_vect[i] = mean(lesions[parcels==i])*100 # get percent of parcel damaged
      pcd_vol[parcels==i] = pcd_vect[i] # put into volume
    }

    pcd_vect=data.frame(Parcel=cfg$node_label,
                        ParcelGroup=cfg$node_group,
                        PercentDamage=pcd_vect)

    if(verbose==T){cat('Finished computing parcel damage measures')}

    # save parcel damage measures

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

    writenii(lesions_or, paste0(pat.path,"/",cfg$pat_id,"_lesion.nii.gz"))
    writenii(pcd_vol, paste0(pd.path,"/",cfg$pat_id,"_",
                             cfg$file_suffix,"_percent_damage.nii.gz"))
    write.csv(pcd_vect,paste0(pd.path,"/",cfg$pat_id,"_",
                              cfg$file_suffix,"_percent_damage.csv"))

  }

}
