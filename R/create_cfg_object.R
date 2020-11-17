#' @title Create cfg structure for analysis
#' @description This function compiles the relevant settings and file paths to be fed into analysis functions.
#' @param pat_id a string specifying the given patient's ID (used as prefix for output files)
#' @param lesion_path a string specifying the path to your lesion mask (pre-registered to MNI template)
#' @param out_path a string specifying the path to your desired output directory; patient and atlas result directories will be created here
#' @param dsi_path a string specifying the path to your local DSI_Studio program
#' @param parcel_path a string specifying the path to a parcellation map (should have identical dimensions to lesion and be in MNI template space);
#' function defaults to the Schaefer-Yeo 100 parcel 7 network parcellation with subcortical structures
#' @param file_suffix a string providing a suffix for results files; the atlas name is recommended (e.g. AAL, Power, Gordon, etc.)
#' @param con_type connectivity type (must be either 'end' or 'pass'): if 'end', connections are defined based on streamline endpoints;
#' if 'pass', they are based on streamline pass-throughs ('end' is more conservative and is recommended)
#' @param sspl_spared_thresh a number between 1 and 100 representing the percent spared threshold for computing shortest structural path lengths
#' (e.g. 100 means that only fully spared regions will be included in SSPL calculation, 1 means that regions with at least 1% spared will be included; default is 50)
#' @param node_label a vector of strings corresponding to node labels (i.e. parcel names); must be the same length (P) as the number of parcels
#' @param node_color a vector of integer values corresponding to e.g. network assignments or partitions (used to color nodes in visualizations)
#' @param parcel_coords a P-by-3 matrix of parcel coordinates used for brain network graphs; if not supplied, they will be estimated from the parcel file
#' @param delta_sspl_thresh a number between 1 and 100 giving the percentile threshold for displaying SSPL increases in figures (will only display SSPL increases above percentile threshold; default = 90)
#' @param parcel_dmg_thresh a number between 1 and 100 giving the parcel damage threshold for displaying SSPL increases in figures (will not display SSPL increases for parcels with % damage greater than or equal to the threshold; default = 100)
#' @param tract_sdc_thresh a number between 1 and 100 giving the tract disconnection threshold for displaying tracts in figures (will not display data for tracts with % disconnection below threshold; default = 5)
#' @param smooth a number representing the sigma (in mm) of the Gaussian kernel for smoothing the track density imaging outputfile (default = 2)
#'
#' @return A list structure to be input into downstream analysis functions.
#' @examples \dontrun{
#'
#' }
#' @export

create_cfg_object=function(pat_id,lesion_path,out_path,
                           dsi_path=NULL,parcel_path=NULL,file_suffix=NULL,
                           con_type="end",sspl_spared_thresh=50,node_label=NULL,
                           node_color=NULL,parcel_coords=NULL,delta_sspl_thresh=90,
                           parcel_dmg_thresh=100,tract_sdc_thresh=5,smooth=2){
  cfg=list()
  cfg$pat_id=pat_id; cfg$lesion_path=lesion_path; cfg$out_path=out_path
  cfg$source_path=system.file("extdata","Tractography_Atlas",package="LQT")

  if(!is.null(dsi_path)){
    if(file.exists(dsi_path)){
      cfg$dsi_path=dsi_path
    }else{
      stop("Specified 'dsi_path' does not exist.")
    }
  }else{
    checkdsi='/Applications/dsi_studio.app/Contents/MacOS/dsi_studio'
    if(file.exists(checkdsi)){
      cfg$dsi_path=checkdsi
    }else{
      stop("Cannot find path to dsi_studio, please specify in function call.")
    }
  }

  if(!is.null(parcel_path)){
    if(file.exists(parcel_path) & !grepl("extdata",parcel_path)){
      cfg$parcel_path=parcel_path
    }else if(file.exists(parcel_path) & grepl("extdata",parcel_path)){
      cfg$parcel_path=parcel_path
      t=read.csv(gsub(".nii.gz",".csv",parcel_path),header=T)
      node_label = t$RegionName
      node_color = t$NetworkID
      parcel_coords = cbind(t$X, t$Y, t$Z)
    }else{
      stop("Specified 'parcel_path' does not exist.")
    }
  }else{
    cfg$parcel_path=system.file("extdata","Schaefer_Yeo_Plus_Subcort",
                                "100Parcels7Networks.nii.gz",package="LQT")
    cat("'parcel_path' not specified; defaulting to Schaefer-Yeo's 100 parcel + 7 network + subcortical structures parcellation.")
    t=read.csv(gsub(".nii.gz",".csv",cfg$parcel_path),header=T)
    node_label = t$RegionName
    node_color = t$NetworkID
    parcel_coords = cbind(t$X, t$Y, t$Z)
  }

  if(!is.null(file_suffix)){
    cfg$file_suffic=file_suffix
  }else if(is.null(file_suffix) & is.null(parcel_path)){
    cfg$file_suffix="Yeo7100"
  }else{
    cfg$file_suffix=tail(strsplit(test,"/|\\.nii|\\.nii\\.gz")[[1]],1)
  }

  if(!con_type%in%c("end","pass")){
    stop("'con_type' must be either 'end' or 'pass'.")
  }

  cfg$con_type=con_type; cfg$sspl_spared_thresh=sspl_spared_thresh
  cfg$node_label=node_label; cfg$node_color=node_color; cfg$smooth=smooth
  cfg$parcel_coords=parcel_coords; cfg$delta_sspl_thresh=delta_sspl_thresh
  cfg$parcel_dmg_thresh=parcel_dmg_thresh; cfg$tract_sdc_thresh=tract_sdc_thresh

  return(cfg)

}

