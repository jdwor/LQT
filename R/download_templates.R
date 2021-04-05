#' @keywords internal

download_templates=function(){
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
