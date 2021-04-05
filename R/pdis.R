#' @keywords internal

pdis = function(x, atlas){
  vars = rownames(atlas)
  pat = apply(x, 1, sum)
  atlas = apply(atlas, 1, sum)

  parc.discon = matrix(pat / atlas,nrow=1)
  colnames(parc.discon)=vars
  return(parc.discon)
}
