#' @title Get tract-based disconnection
#' @description This function computes tract-based disconnection measures using an MNI-registered lesion and the tract segmentations
#' obtained from the curated HCP-842 tractography atlas as described in Yeh et al., (2018 - NeuroImage).
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

get_tract_discon<-function(cfg,saveout=TRUE){
  tract_path = paste0(cfg$source_path,"/All_Tracts")
  my_tracts = list.files(tract_path, pattern="\\.trk\\.gz$")

  pat.path=paste0(cfg$out_path,"/",cfg$pat_id)
  if(!dir.exists(cfg$out_path)){
    dir.create(cfg$out_path)
  }
  if(!dir.exists(pat.path)){
    dir.create(pat.path)
  }

  td.path=paste0(pat.path,"/Tract_Disconnection")
  if(!dir.exists(td.path)){
    dir.create(td.path)
  }

  for(i in 1:length(my_tracts)){
    # output file name
    tract_name = substr(my_tracts[i],1,nchar(my_tracts[i])-7)
    out_file = paste0(td.path,"/",my_tracts[i])

    # compute tract disconnection
    system(paste0("! ",cfg$dsi_path,' --action=ana --source=',cfg$source_path,"/HCP842_1mm.fib.gz",
                  " --tract=",cfg$source_path,"/All_Tracts/",my_tracts[i]," --roi=",cfg$lesion_path,
                  " --output=",out_file," --export=stat"))

    # Create empty text file if it does not exist (necessary because some versions of DSI_Studio don't output anything for null results)
    if(!file.exists(paste0(td.path,"/",tract_name,'.trk.gz.stat.txt'))){
      fileConn<-file(paste0(td.path,"/",tract_name,'.trk.gz.stat.txt'))
      writeLines("number of tracts\t0", fileConn)
      close(fileConn)
    }
  }

  to_remove=list.files(td.path, pattern="\\.trk\\.gz$")
  file.remove(paste0(td.path,"/",to.remove))

  my_stats=list.files(td.path, pattern="\\.txt$")

  tract_name=rep(NA,length(my_stats))
  tract_discon=rep(NA,length(my_stats))
  for(i in 1:length(my_stats)){
    fid=read.table(paste0(td.path,"/",my_stats[i]),sep="\t")

    tract_name[i] = substr(my_stats[i],1,nchar(my_tracts[i])-16)
    tract_discon[i] = fid$V2[1]
  }

  to_remove=list.files(td.path, pattern="\\.txt$")
  file.remove(paste0(td.path,"/",to.remove))

  tract_names = readMat(paste0(tract_path,"/","tract_names.mat"))
  tract_names = as.character(unlist(tract_names))
  tract_counts = readMat(paste0(tract_path,"/","tract_counts.mat"))
  tract_counts = as.numeric(unlist(tract_counts))

  if(identical(tract_names,tract_name)){
    tract_discon = tract_discon/tract_counts
  }else{
    tc_index = match(tract_name,tract_names)
    tract_discon = 100*(tract_discon/tract_counts[tc_index])
  }

  cat("Finished computing patient disconnection measures.")

  output = data.frame(Tract = tract_name, Discon = tract_discon)

  if(saveout==T){
    write.csv(output,paste0(td.path,cfg$pat_id,"_percent_discon_tracts.csv"))
  }

  return(output)
}
