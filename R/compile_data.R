#' @title Compile datasets
#' @description This function takes completed LQT files and compiles analysis-ready datasets
#' @param out_path a string specifying the output directory to which LQT files were saved.
#'
#' @importFrom neurobase readnii writenii
#' @importFrom utils read.csv
#'
#' @return An .RData file with the suffix .connectivity.RData. This contains the structural disconnection matrix (connectivity);
#' an .RData file with the suffix .network_measures.RData, which contains various graph measures for the SC matrix;
#' an .RData file with the suffix _percent_parcel_mats.RData. This file contains a disconnection adjacency matrix (pct_sdc_matrix) and a spared connection adjacency matrix (pct_spared_sc_matrix);
#'
#' @export

compile_data<-function(out_path){
  subjects=list.files(out_path)
  at.path=file.path(out_path,"Atlas")
  pat_ids=file.path(out_path,list.files(out_path))


}
