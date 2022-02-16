#' @title Create configuration (cfg) structure for analysis
#' @description This function compiles the relevant settings and file paths to be fed into analysis functions.
#' @param pat_ids a vector of strings specifying patients' IDs (used as directories for output files)
#' @param lesion_paths a vector of strings specifying the path to patients' lesion masks (pre-registered to MNI template, and corresponding to pat_ids)
#' @param out_path a string specifying the path to your desired output directory; patient and atlas result directories will be created here
#' @param dsi_path a string specifying the path to your local DSI_Studio program
#' @param parcel_path a string specifying the path to a parcellation map (should have identical dimensions to lesion and be in MNI template space);
#' function defaults to the Schaefer-Yeo 100 parcel 7 network parcellation with subcortical structures
#' @param file_suffix a string providing a suffix for results files; the atlas name is recommended (e.g. AAL, Power, Gordon, etc.)
#' @param con_type connectivity type (must be either 'end' or 'pass'): if 'end', connections are defined based on streamline endpoints;
#' if 'pass', they are based on streamline pass-throughs ('end' is more conservative and is recommended)
#' @param sspl_spared_thresh a number between 1 and 100 representing the percent spared threshold for computing shortest structural path lengths
#' (e.g. 100 means that only fully spared regions will be included in SSPL calculation, 1 means that regions with at least 1 percent spared will be included; default is 50)
#' @param node_label a vector of strings corresponding to node labels (i.e. parcel names); must be the same length (P) as the number of parcels
#' @param node_group a vector of strings integer values corresponding to e.g. network assignments or partitions (used to color nodes in visualizations)
#' @param parcel_coords a P-by-3 matrix of parcel coordinates used for brain network graphs; if not supplied, they will be estimated from the parcel file
#' @param delta_sspl_thresh a number between 1 and 100 giving the percentile threshold for displaying SSPL increases in figures (will only display SSPL increases above percentile threshold; default = 90)
#' @param parcel_dmg_thresh a number between 1 and 100 giving the parcel damage threshold for displaying SSPL increases in figures (will not display SSPL increases for parcels with percent damage greater than or equal to the threshold; default = 100)
#' @param tract_sdc_thresh a number between 1 and 100 giving the tract disconnection threshold for displaying tracts in figures (will not display data for tracts with percent disconnection below threshold; default = 5)
#' @param smooth a number representing the smoothing kernel FWHM (in voxel units) for smoothing the track density imaging output file (default = 2)
#'
#' @importFrom utils read.csv tail download.file unzip
#' @importFrom purrr map list_modify
#' @importFrom dplyr %>%
#' @importFrom neurobase check_mask
#'
#' @return A list structure to be input into downstream analysis functions.
#'
#' @export

create_cfg_object=function(pat_ids,lesion_paths,out_path,
                           dsi_path=NULL,parcel_path=NULL,file_suffix=NULL,
                           con_type="end",sspl_spared_thresh=50,node_label=NULL,
                           node_group=NULL,parcel_coords=NULL,delta_sspl_thresh=90,
                           parcel_dmg_thresh=100,tract_sdc_thresh=5,smooth=2){
  cfg=list()
  cfg$pat_id=pat_ids; cfg$lesion_path=lesion_paths
  if(length(pat_ids)!=length(lesion_paths)){
    stop("Number of patient ids and lesion paths must be the same.")
  }else{
    num_subs=length(pat_ids)
  }
  cfg$out_path=rep(out_path,num_subs)

  if(!file.exists(system.file("extdata","Tractography_Atlas",package="LQT"))){
    download_templates()
  }

  cfg$source_path=rep(system.file("extdata","Tractography_Atlas",package="LQT"),num_subs)

  if(!is.null(dsi_path)){
    if(file.exists(dsi_path)){
      cfg$dsi_path=rep(dsi_path,num_subs)
    }else{
      stop("Specified path to dsi_studio ('dsi_path') does not exist.\n
           If it is already installed, please specify the correct path in the function call.\n
           If not, please download the software at http://dsi-studio.labsolver.org/dsi-studio-download.
           For easiest operability of LQT, please download the 08/30/21 version of DSI_studio.")
    }
  }else{
    if(.Platform$OS.type=="unix"){
      checkdsi=system.file("extdata","DSI_studio/dsi_studio.app/Contents/MacOS/dsi_studio",
                           package="LQT")
    }else if(.Platform$OS.type=="windows"){
      checkdsi=system.file("extdata","DSI_studio/dsi_studio_64/dsi_studio.exe",
                           package="LQT")
    }else{
      checkdsi="dsi_path not found"
    }
    if(file.exists(checkdsi)){
      cfg$dsi_path=rep(checkdsi,num_subs)
    }else{
      stop("Cannot find path to dsi_studio.\n
           If it is already installed, please specify the correct path in the function call.\n
           If not, please download the software at http://dsi-studio.labsolver.org/dsi-studio-download.
           For easiest operability of LQT, please download the 08/30/21 version of DSI_studio.")
    }
  }

  if(!is.null(parcel_path)){
    if(file.exists(parcel_path) & !grepl("extdata",parcel_path)){
      cfg$parcel_path=rep(parcel_path,num_subs)
      if(is.null(node_label)){
        warning("No node labels specified; nodes will be labeled numerically.")
      }
      if(is.null(node_group)){
        warning("No node groups specified.")
      }
      if(is.null(parcel_coords)){
        warning("No parcel coordinates specified; network plots will not be possible")
      }
      node_label = t$RegionName
      node_group = t$NetworkID
      parcel_coords = cbind(t$X, t$Y, t$Z)
    }else if(file.exists(parcel_path) & grepl("extdata",parcel_path)){
      cfg$parcel_path=rep(parcel_path,num_subs)
      t=read.csv(gsub(".nii.gz",".csv",parcel_path[1]),header=T)
      node_label = t$RegionName
      node_group = t$NetworkID
      parcel_coords = cbind(t$X, t$Y, t$Z)
    }else if(parcel_path==""){
      stop("\n\n'parcel_path' is empty. You might have used a system.file call before templates were downloaded. If so, re-run the parcel_path call and try again.")
    }else{
      stop("Specified 'parcel_path' does not exist.")
    }
  }else{
    cfg$parcel_path=rep(system.file("extdata","Schaefer_Yeo_Plus_Subcort",
                        "100Parcels7Networks.nii.gz",package="LQT"),num_subs)
    cat("'parcel_path' not specified; defaulting to Schaefer-Yeo's 100 parcel + 7 network + subcortical structures parcellation.")
    t=read.csv(gsub(".nii.gz",".csv",cfg$parcel_path[1]),header=T)
    node_label = t$RegionName
    node_group = t$NetworkID
    parcel_coords = cbind(t$X, t$Y, t$Z)
  }

  if(!is.null(file_suffix)){
    cfg$file_suffix=rep(file_suffix,num_subs)
  }else if(is.null(file_suffix) & is.null(parcel_path)){
    cfg$file_suffix=rep("Yeo.100.7",num_subs)
  }else if(is.null(file_suffix) & grepl("Schaefer_Yeo",parcel_path)){
    parc_nums=tail(strsplit(parcel_path,"/|\\.nii|\\.nii\\.gz")[[1]],1)
    parc_nums=parc_nums %>% gsub("Parcels","\\.",.) %>%
      gsub("Networks","",.) %>% paste0("Yeo.",.)
    cfg$file_suffix=rep(parc_nums,num_subs)
  }else{
    cfg$file_suffix=rep(tail(strsplit(parcel_path,"/|\\.nii|\\.nii\\.gz")[[1]],1),num_subs)
  }

  if(!con_type%in%c("end","pass")){
    stop("'con_type' must be either 'end' or 'pass'.")
  }

  cfg$con_type=rep(con_type,num_subs)
  cfg$sspl_spared_thresh=rep(sspl_spared_thresh,num_subs)
  cfg$smooth=rep(smooth,num_subs)
  cfg$delta_sspl_thresh=rep(delta_sspl_thresh,num_subs)
  cfg$parcel_dmg_thresh=rep(parcel_dmg_thresh,num_subs)
  cfg$tract_sdc_thresh=rep(tract_sdc_thresh,num_subs)

  if(num_subs>1){
    cfg = lapply(1:num_subs,function(x, cfg){
      cfg %>% map(x) %>% list_modify(node_label = node_label,
                                     node_group = node_group,
                                     parcel_coords = parcel_coords)
      },cfg)
  }else{
    cfg$node_label=node_label
    cfg$node_group=node_group
    cfg$parcel_coords=parcel_coords
  }

  lesion1=readnii(lesion_paths[1])
  is_mask = check_mask(lesion1)
  if(!is_mask){
    lesion_paths_bin=gsub("\\.nii","\\_bin\\.nii",lesion_paths)
    warning("Lesion masks do not appear to be binary!",immediate. = T)
    thresh=readline(prompt="Enter a numeric threshold to binarize the masks, or enter 'c' to cancel: ")
    if(thresh=='c'){
      stop("Function cancelled, no cfg object created.\nPlease create binary lesion masks and try again.")
    }else if(!is.na(as.numeric(thresh))){
      thresh=as.numeric(thresh)
      for(i in 1:length(lesion_paths)){
        lesion=readnii(lesion_paths[i])
        lesion[lesion<=thresh]=0
        lesion[lesion>thresh]=1
        writenii(lesion,lesion_paths_bin[i])
      }
      cfg$lesion_path=lesion_paths_bin
    }else{
      stop("Must enter either a numeric threshold or 'c'.")
    }
  }

  return(cfg)

}

