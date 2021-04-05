#' @keywords internal

p2pdis = function(x, atlas){
  pat = x
  pat[upper.tri(pat)]=NA
  pat = melt(pat, na.rm=T)
  pat$Var = paste0(pat$Var2, "_to_",pat$Var1)

  atlas[upper.tri(atlas)]=NA
  atlas = melt(atlas, na.rm=T)
  atlas$Var = paste0(atlas$Var2, "_to_",atlas$Var1)

  pat$value = pat$value / atlas$value
  parc2parc.discon = matrix(pat$value,nrow=1)
  colnames(parc2parc.discon)=pat$Var
  return(parc2parc.discon)
}
