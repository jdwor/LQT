#' @keywords internal

ndis = function(x, atlas, groups){

  vars = rownames(atlas)
  pat = apply(x, 1, sum)
  atlas = apply(atlas, 1, sum)

  pat=rowsum(pat, group=groups)
  atlas=rowsum(atlas, group=groups)

  net.discon = matrix(pat / atlas,nrow=1)
  colnames(net.discon)=rownames(pat)
  return(net.discon)
}
