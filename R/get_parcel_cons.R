#' @title Get parcel disconnection measures/maps
#' @description This is a wrapper function that calls a series of functions to create various parcel disconnection measures.
#' It wraps the 'get_parcel_atlas', 'get_atlas_sspl', 'get_parcel_discon', and 'get_patient_sspl' functions.
#' @param cfg a pre-made cfg structure (as list object).
#' @param cores an integer value that indicates how many parallel cores the function should be run on.
#'
#' @export

get_parcel_cons<-function(cfg, cores=1){

  if(is.null(cfg$pat_id)){
    acfg=cfg[[1]]
    at.path=paste0(acfg$out_path,"/Atlas")
  }else{
    acfg=cfg
    at.path=paste0(acfg$out_path,"/Atlas")
  }

  if(!dir.exists(at.path)){
    # get parcel SC for atlas
    get_parcel_atlas(acfg)
    get_atlas_sspl(acfg)
  }else{
    # check whether atlas connectivity matrix exists
    con.file=list.files(at.path,pattern=paste0(acfg$file_suffix,"_",
                                               "connectivity\\.RData$"))
    if(length(con.file)==0){
      get_parcel_atlas(acfg) # if not, create it
    }

    # check whether atlas sspl file exists
    sspl.file=list.files(at.path,pattern=paste0(acfg$file_suffix,"_",
                                                "SSPL_matrix\\.RData$"))
    if(length(sspl.file)==0){
      get_atlas_sspl(acfg) # if not, create it
    }
  }

  # get parcel SDC for patient
  get_parcel_discon(cfg, cores=cores)

  # get parcel SSPL and delta SSPL for patient
  get_patient_sspl(cfg, cores=cores)
}


