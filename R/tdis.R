#' @keywords internal

tdis = function(x){
  td.file = list.files(file.path(x, "Tract_Disconnection"))
  td.file = td.file[grepl(".csv",td.file)]
  td.path = file.path(x, "Tract_Disconnection", td.file)
  tract_discon = read.csv(td.path)
  tract.discon = matrix(tract_discon$Discon,nrow=1)
  colnames(tract.discon)=tract_discon$Tract
  return(tract.discon)
}
