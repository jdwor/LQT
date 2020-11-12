#' @title Get parcel disconnection measures/maps
#' @description This is a wrapper function that calls a series of functions to create various parcel disconnection measures.
#' It wraps the 'get_parcel_atlas', 'get_atlas_sspl', 'get_parcel_discon', and 'get_patient_sspl' functions.
#' @param cfg a pre-made cfg structure (as list object).
#'
#' @examples \dontrun{
#'
#' }
#' @export

get_parcel_cons<-function(cfg){
  at.path=paste0(cfg$out_path,"/Atlas")
  if(!dir.exists(at.path)){
    # get parcel SC for atlas
    get_parcel_atlas(cfg)
    get_atlas_sspl(cfg)
  }else{
    # check whether atlas connectivity matrix exists
    con.file=list.files(at.path,pattern="connectivity\\.RData$")
    if(length(con.file)==0){
      get_parcel_atlas(cfg) # if not, create it
    }

    # check whether atlas sspl file exists
    sspl.file=list.files(at.path,pattern="SSPL_matrix\\.RData$")
    if(length(sspl.file)==0){
      get_atlas_sspl(cfg) # if not, create it
    }
  }

  # get parcel SDC for patient
  get_parcel_discon(cfg)

  # get parcel SSPL and delta SSPL for patient
  get_patient_sspl(cfg)
}


