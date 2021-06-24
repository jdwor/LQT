#' @title Get tract-based disconnection
#' @description This function computes tract-based disconnection measures using an MNI-registered lesion and the tract segmentations
#' obtained from the curated HCP-842 tractography atlas as described in Yeh et al., (2018 - NeuroImage).
#' @param cfg a pre-made cfg structure (as list object).
#' @param cores an integer value that indicates how many parallel cores the function should be run on.
#'
#' @importFrom neurobase readnii writenii
#' @importFrom R.matlab readMat
#' @importFrom utils read.table write.csv
#' @importFrom pbmcapply pbmclapply
#'
#' @return A .csv file giving the percent of streamlines in each tract that were disconnected by the lesions (column "Discon")
#' and the associated tract names (column "Tract").
#'
#' @export

get_tract_discon<-function(cfg, cores=1){
  if(is.null(cfg$pat_id)){

    cat('Computing tract disconnection.\n')
    out = pbmclapply(cfg, get_tract_discon, mc.cores=cores)
    cat("Finished computing tract disconnection measures.")

  }else{

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

    cat('Computing tract disconnection.\n')

    num_tracts=length(my_tracts)
    for(i in 1:num_tracts){
      # output file name
      tract_name = substr(my_tracts[i],1,nchar(my_tracts[i])-7)
      out_file = paste0(td.path,"/",my_tracts[i])

      # compute tract disconnection
      out=tryCatch({
        suppressMessages(system(paste0("! ",cfg$dsi_path," --action=ana --source=",cfg$source_path,"/HCP842_1mm.fib.gz",
                                       " --tract=",cfg$source_path,"/All_Tracts/",my_tracts[i]," --roi=",cfg$lesion_path,
                                       " --output=",out_file," --export=stat"),intern=T))
      },error=function(e){
        suppressMessages(shell(paste0(cfg$dsi_path," --action=ana --source=",cfg$source_path,"/HCP842_1mm.fib.gz",
                                       " --tract=",cfg$source_path,"/All_Tracts/",my_tracts[i]," --roi=",cfg$lesion_path,
                                       " --output=",out_file," --export=stat"),intern=T))
      })

      # Create empty text file if it does not exist (necessary because some versions of DSI_Studio don't output anything for null results)
      if(!file.exists(paste0(td.path,"/",tract_name,'.trk.gz.stat.txt'))){
        fileConn<-file(paste0(td.path,"/",tract_name,'.trk.gz.stat.txt'))
        writeLines("number of tracts\t0", fileConn)
        close(fileConn)
      }else if(file.info(paste0(td.path,"/",tract_name,'.trk.gz.stat.txt'))$size==0){
        fileConn<-file(paste0(td.path,"/",tract_name,'.trk.gz.stat.txt'))
        writeLines("number of tracts\t0", fileConn)
        close(fileConn)
      }

      cat('Progress:',i,'of',num_tracts,'tracts evaluated\r')
    }

    to_remove=list.files(td.path, pattern="\\.trk\\.gz$")
    if(length(to_remove)>0){file.remove(paste0(td.path,"/",to_remove))}

    my_stats=list.files(td.path, pattern="\\.txt$")

    tract_name=rep(NA,length(my_stats))
    tract_discon=rep(NA,length(my_stats))
    for(i in 1:length(my_stats)){
      fid=read.table(paste0(td.path,"/",my_stats[i]),sep="\t")

      tract_name[i] = substr(my_stats[i],1,nchar(my_stats[i])-16)
      tract_discon[i] = fid$V2[1]
    }

    to_remove=list.files(td.path, pattern="\\.txt$")
    file.remove(paste0(td.path,"/",to_remove))

    load(paste0(tract_path,"/","tract_info.RData"))
    if(identical(tract_names,tract_name)){
      tract_discon = tract_discon/tract_counts
    }else{
      tc_index = match(tract_name,tract_names)
      tract_discon = 100*(tract_discon/tract_counts[tc_index])
      tract_pathways = tract_pathways[tc_index]
    }

    output = data.frame(Tract = tract_name, Discon = tract_discon,
                        Pathway = tract_pathways)
    write.csv(output,paste0(td.path,"/",cfg$pat_id,"_percent_discon_tracts.csv"))

    cat("Finished computing tract disconnection measures.")

  }
}
