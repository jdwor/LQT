#' @title Compile datasets
#' @description This function takes completed LQT files and compiles analysis-ready datasets
#' @param cfg a pre-made cfg structure (as list object).
#' @param out_paths a string specifying the output directory to which LQT files were saved,
#' or a vector of strings specifying individual patients' output directories.
#'
#' @importFrom neurobase readnii writenii
#' @importFrom utils read.csv
#' @importFrom reshape2 melt
#'
#' @return An .RData file with the suffix .connectivity.RData. This contains the structural disconnection matrix (connectivity);
#' an .RData file with the suffix .network_measures.RData, which contains various graph measures for the SC matrix;
#' an .RData file with the suffix _percent_parcel_mats.RData. This file contains a disconnection adjacency matrix (pct_sdc_matrix) and a spared connection adjacency matrix (pct_spared_sc_matrix);
#'
#' @export

compile_data<-function(cfg=NULL, out_paths=NULL, parallel=F, cores=2){

  if(!is.null(cfg) & !is.list(cfg)){
    out_paths = cfg; cfg = NULL
  }
  if(is.null(cfg) & is.null(out_paths)){
    stop("Either 'cfg' or 'out_paths' must be specified.")
  }else if(!is.null(cfg)){
    pat_ids=unlist(lapply(cfg,`[[`,1))
    out_path=cfg[[1]]$out_path
    at.path=file.path(out_path,"Atlas")
    pat.paths=file.path(out_path,pat_ids)
  }else{
    if(length(out_paths)==1){
      pat_ids=list.dirs(out_path,recursive=F)
      pat_ids=pat_ids[pat_ids!="Atlas"]
      at.path=file.path(out_path,"Atlas")
      pat.paths=file.path(out_path,pat_ids)
    }else{
      which.atlas=grepl("/Atlas",out_paths)
      if(sum(which.atlas)==0){stop("Path to 'Atlas' folder must be included in 'out_paths.")}
      at.path=out_paths[which.atlas]
      pat.paths=out_paths[!which.atlas]
    }
  }

  ID = pat_ids
  at.file = list.files(at.path)
  at.file = at.file[grepl("connectivity.RData",at.file)]
  load(file.path(at.path,at.file))
  atlas = connectivity; rm(connectivity)

  parc.damage = ifelse(parallel=T,lapply(pat.paths,pdam),
                       mclapply(pat.paths,pdam,cores=cores))
  parc.damage = do.call(rbind, parc.damage)
  parc.damage = as.data.frame(parc.damage)
  parc.damage = cbind(ID, parc.damage)

  tract.discon = ifelse(parallel=T,lapply(pat.paths,tdis),
                        mclapply(pat.paths,tdis,cores=cores))
  tract.discon = do.call(rbind, tract.discon)
  tract.discon = as.data.frame(tract.discon)
  tract.discon = cbind(ID, tract.discon)

  patnets = ifelse(parallel=T,lapply(pat.paths,p2pnet),
                   mclapply(pat.paths,p2pnet,cores=cores))

  net.discon = ifelse(parallel=T,lapply(patnets,ndis,node_group),
                      mclapply(patnets,ndis,node_group,cores=cores))
  net.discon = do.call(rbind, net.discon)
  net.discon = as.data.frame(net.discon)
  net.discon = cbind(ID, net.discon)

  parc.discon = ifelse(parallel=T,lapply(patnets,pdis),
                       mclapply(patnets,pdis,cores=cores))
  parc.discon = do.call(rbind, parc.discon)
  parc.discon = as.data.frame(parc.discon)
  parc.discon = cbind(ID, parc.discon)

  net2net.discon = ifelse(parallel=T,lapply(patnets,n2ndis,node_group),
                          mclapply(patnets,n2ndis,node_group,cores=cores))
  net2net.discon = do.call(rbind, net2net.discon)
  net2net.discon = as.data.frame(net2net.discon)
  net2net.discon = cbind(ID, net2net.discon)

  parc2parc.discon = ifelse(parallel=T,lapply(patnets,p2pdis),
                            mclapply(patnets,p2pdis,cores=cores))
  parc2parc.discon = do.call(rbind, parc2parc.discon)
  parc2parc.discon = as.data.frame(parc2parc.discon)
  parc2parc.discon = cbind(ID, parc2parc.discon)

  output = list(parc.damage = parc.damage, tract.discon = tract.discon,
                net.discon = net.discon, parc.discon = parc.discon,
                net2net.discon = net2net.discon,
                parc2parc.discon = parc2parc.discon)
  return(output)

}
