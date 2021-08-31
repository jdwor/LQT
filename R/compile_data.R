#' @title Compile datasets
#' @description This function takes completed LQT files and compiles analysis-ready datasets
#' @param cfg a pre-made cfg structure (as list object).
#' @param cores an integer value that indicates how many parallel cores the function should be run on.
#'
#' @importFrom neurobase readnii writenii
#' @importFrom utils read.csv
#' @importFrom reshape2 melt
#' @importFrom parallel mclapply
#'
#' @return A list object containing several analysis-ready datasets, in which rows represent patients and columns represent lesion metrics:
#' "parc.damage", which gives percent regional parcel damage across subjects;
#' "tract.discon", which gives percent tract disconnection across subjects;
#' "net.discon", which gives the percent disconnection for each parcel group or network across its edges;
#' "parc.discon", which gives the percent disconnection for each individual parcel across its edges;
#' "net2net.discon", which gives the percent disconnection of each pairwise network-to-network edge;
#' "parc2parc.discon", which gives the percent disconnection of each pairwise parcel-to-parcel edge.
#'
#' @export

compile_data<-function(cfg, cores=1){

  pat_ids=unlist(lapply(cfg,`[[`,1))
  out_path=cfg[[1]]$out_path
  suffix=cfg[[1]]$file_suffix
  at.path=file.path(out_path,"Atlas")
  pat.paths=file.path(out_path,pat_ids)

  ID = pat_ids
  load(paste0(at.path,'/atlas_',suffix,'_connectivity.RData'))
  atlas = connectivity; rm(connectivity)

  parc.damage = mclapply(pat.paths,pdam,suffix,mc.cores=cores)
  parc.damage = do.call(rbind, parc.damage)
  parc.damage = as.data.frame(parc.damage)
  parc.damage = cbind(ID, parc.damage)

  tract.discon = mclapply(pat.paths,tdis,mc.cores=cores)
  tract.discon = do.call(rbind, tract.discon)
  tract.discon = as.data.frame(tract.discon)
  tract.discon = cbind(ID, tract.discon)

  patnets = mclapply(pat.paths,p2pnet,suffix,mc.cores=cores)

  net.discon = mclapply(patnets,ndis,atlas,node_group,mc.cores=cores)
  net.discon = do.call(rbind, net.discon)
  net.discon = as.data.frame(net.discon)
  net.discon = cbind(ID, net.discon)

  parc.discon = mclapply(patnets,pdis,atlas,mc.cores=cores)
  parc.discon = do.call(rbind, parc.discon)
  parc.discon = as.data.frame(parc.discon)
  parc.discon = cbind(ID, parc.discon)

  net2net.discon = mclapply(patnets,n2ndis,atlas,node_group,mc.cores=cores)
  net2net.discon = do.call(rbind, net2net.discon)
  net2net.discon = as.data.frame(net2net.discon)
  net2net.discon = cbind(ID, net2net.discon)

  parc2parc.discon = mclapply(patnets,p2pdis,atlas,mc.cores=cores)
  parc2parc.discon = do.call(rbind, parc2parc.discon)
  allnas = apply(parc2parc.discon, 2, function(x) mean(is.na(x))==1)
  parc2parc.discon = as.data.frame(parc2parc.discon[,!allnas])
  parc2parc.discon = cbind(ID, parc2parc.discon)

  output = list(parc.damage = parc.damage, tract.discon = tract.discon,
                net.discon = net.discon, parc.discon = parc.discon,
                net2net.discon = net2net.discon,
                parc2parc.discon = parc2parc.discon)
  return(output)

}
