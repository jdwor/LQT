#' @keywords internal

download_templates=function(){
  extdata=system.file("extdata",package="LQT")
  cat("Template files have not yet been installed. Downloading them now from Figshare.\nThis may take a few minutes, but will only happen once...\n")
  download.file("https://ndownloader.figshare.com/files/27368315?private_link=2d830ec228a1c4bdf8aa",
                destfile=paste0(extdata,"/HCP842_QA.nii.gz"),mode="wb")
  download.file("https://ndownloader.figshare.com/files/27368318?private_link=2d830ec228a1c4bdf8aa",
                destfile=paste0(extdata,"/MNI152_T1_1mm.nii.gz"),mode="wb")
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
  if(.Platform$OS.type=="unix"){
    download.file("https://figshare.com/ndownloader/files/31708106?private_link=c67b652f2056867a5648",
                  destfile=paste0(extdata,"/DSI_studio.zip"),mode="wb")
    unzip(zipfile = paste0(extdata,"/DSI_studio.zip"),
          exdir = paste0(extdata,"/DSI_studio"))
    file.remove(paste0(extdata,"/DSI_studio.zip"))
  }else if(.Platform$OS.type=="windows"){
    download.file("https://figshare.com/ndownloader/files/31708160?private_link=14b4549b8479c88a8343",
                  destfile=paste0(extdata,"/DSI_studio.zip"),mode="wb")
    unzip(zipfile = paste0(extdata,"/DSI_studio.zip"),
          exdir = paste0(extdata,"/DSI_studio"))
    file.remove(paste0(extdata,"/DSI_studio.zip"))
  }
}
