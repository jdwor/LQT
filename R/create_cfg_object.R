#' @title Create cfg structure for analysis
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
#' @param smooth a number representing the sigma (in mm) of the Gaussian kernel for smoothing the track density imaging output file (default = 2)
#'
#' @importFrom utils read.csv tail download.file unzip
#' @importFrom purrr map list_modify
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
    extdata=system.file("extdata",package="LQT")
    cat("Template files have not yet been installed. Downloading them now from Figshare.\nThis may take a few minutes, but will only happen once...\n")
    download.file("https://ndownloader.figshare.com/files/27368315?private_link=2d830ec228a1c4bdf8aa",
                  destfile=paste0(extdata,"/HCP842_QA.nii.gz"))
    download.file("https://ndownloader.figshare.com/files/27368318?private_link=2d830ec228a1c4bdf8aa",
                  destfile=paste0(extdata,"/MNI152_T1_1mm.nii.gz"))
    download.file("https://ndownloader.figshare.com/articles/14343335?private_link=66e65823e263a46c8237",
                  destfile=paste0(extdata,"/Schaefer_Yeo_Plus_Subcort.zip"),mode="wb")
    unzip(zipfile = paste0(extdata,"/Schaefer_Yeo_Plus_Subcort.zip"),
          exdir = paste0(extdata,"/Schaefer_Yeo_Plus_Subcort"))
    file.remove(paste0(extdata,"/Schaefer_Yeo_Plus_Subcort.zip"))
    download.file("https://ndownloader.figshare.com/articles/14343344?private_link=796fc73e7fc4fa6a1bca",
                  destfile=paste0(extdata,"/Other_Atlases.zip"),mode="wb")
    unzip(zipfile = paste0(extdata,"/Other_Atlases.zip"),
          exdir = paste0(extdata,"/Other_Atlases"))
    file.remove(paste0(extdata,"/Other_Atlases.zip"))
    download.file("https://ndownloader.figshare.com/articles/14342426?private_link=4be1178860a7d8dad555",
                  destfile=paste0(extdata,"/Tractography_Atlas.zip"),mode="wb")
    unzip(zipfile = paste0(extdata,"/Tractography_Atlas.zip"),
          exdir = paste0(extdata,"/Tractography_Atlas"))
    file.remove(paste0(extdata,"/Tractography_Atlas.zip"))
    download.file("https://ndownloader.figshare.com/articles/14342450?private_link=83a8d620899ed9b198d3",
                  destfile=paste0(extdata,"/Tractography_Atlas/All_Tracts.zip"),mode="wb")
    unzip(zipfile = paste0(extdata,"/Tractography_Atlas/All_Tracts.zip"),
          exdir = paste0(extdata,"/Tractography_Atlas/All_Tracts"))
    file.remove(paste0(extdata,"/Tractography_Atlas/All_Tracts.zip"))
  }

  cfg$source_path=rep(system.file("extdata","Tractography_Atlas",package="LQT"),num_subs)

  if(!is.null(dsi_path)){
    if(file.exists(dsi_path)){
      cfg$dsi_path=rep(dsi_path,num_subs)
    }else{
      stop("Specified path to dsi_studio ('dsi_path') does not exist.\n
           If it is already installed, please specify the correct path in the function call.\n
           If not, please download the software at http://dsi-studio.labsolver.org/dsi-studio-download.")
    }
  }else{
    checkdsi='/Applications/dsi_studio.app/Contents/MacOS/dsi_studio'
    if(file.exists(checkdsi)){
      cfg$dsi_path=rep(checkdsi,num_subs)
    }else{
      stop("Cannot find path to dsi_studio.\n
           If it is already installed, please specify the correct path in the function call.\n
           If not, please download the software at http://dsi-studio.labsolver.org/dsi-studio-download.")
    }
  }

  if(!is.null(parcel_path)){
    if(file.exists(parcel_path) & !grepl("extdata",parcel_path)){
      cfg$parcel_path=rep(parcel_path,num_subs)
    }else if(file.exists(parcel_path) & grepl("extdata",parcel_path)){
      cfg$parcel_path=rep(parcel_path,num_subs)
      t=read.csv(gsub(".nii.gz",".csv",parcel_path[1]),header=T)
      node_label = t$RegionName
      node_group = t$NetworkID
      parcel_coords = cbind(t$X, t$Y, t$Z)
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
    cfg$file_suffix=rep("Yeo7100",num_subs)
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

  return(cfg)

}

