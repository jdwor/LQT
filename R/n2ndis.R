#' @keywords internal

n2ndis = function(x, atlas, groups){
  pat = x; rownames(pat)=groups; colnames(pat)=groups
  pat = pat[order(groups),order(groups)]
  pat[upper.tri(pat)]=NA
  pat = melt(pat, na.rm=T)
  pat$Var = paste0(pat$Var2, "_to_",pat$Var1)

  rownames(atlas)=groups; colnames(atlas)=groups
  atlas = atlas[order(groups),order(groups)]
  atlas[upper.tri(atlas)]=NA
  atlas = melt(atlas, na.rm=T)
  atlas$Var = paste0(atlas$Var2, "_to_",atlas$Var1)

  pat = rowsum(pat$value, group=pat$Var)
  atlas = rowsum(atlas$value, group=atlas$Var)
  discon = pat / atlas; discon[is.na(discon)]=0

  net2net.discon = matrix(discon ,nrow=1)
  colnames(net2net.discon)=rownames(pat)
  return(net2net.discon)
}
