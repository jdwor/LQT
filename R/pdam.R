#' @keywords internal

pdam = function(x, suffix){
  pd.file = list.files(file.path(x, "Parcel_Damage"))
  pd.file = pd.file[grepl(paste0(suffix,"_percent_damage\\.csv"),pd.file)]
  pd.path = file.path(x, "Parcel_Damage", pd.file)
  parcel_damage = read.csv(pd.path)
  parcel.damage = matrix(parcel_damage$PercentDamage,nrow=1)
  colnames(parcel.damage)=parcel_damage$Parcel
  return(parcel.damage)
}
