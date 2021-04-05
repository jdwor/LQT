#' @keywords internal

p2pnet = function(x){
  pd.file = list.files(file.path(x, "Parcel_Disconnection"))
  pd.file = pd.file[grepl("connectivity.RData",pd.file)]
  pd.path = file.path(x, "Parcel_Disconnection", pd.file)
  load(pd.path)
  return(connectivity)
}
