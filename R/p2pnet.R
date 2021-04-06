#' @keywords internal

p2pnet = function(x, suffix){
  pd.file = list.files(file.path(x, "Parcel_Disconnection"))
  pd.file = pd.file[grepl(paste0(suffix,"_disconnectivity.RData"),pd.file)]
  pd.path = file.path(x, "Parcel_Disconnection", pd.file)
  load(pd.path)
  dmat = disconnectivity
  return(dmat)
}
