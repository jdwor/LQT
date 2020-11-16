#' @keywords internal
#' @importFrom neurobase readnii writenii

util_get_coords=function(cfg){
  # load parcels
  my_img = readnii(cfg$parcel_path)

  # get unique values
  unique_v = unique(as.vector(my_img));
  unique_v = unique_v[unique_v > 0]

  # loop through unique values and get coordinates in MNI
  get.coords=function(x,unique_v,my_img){
    my_coords = which(my_img==unique_v[x],arr.ind=T)
    av_coords = apply(my_coords,2,mean)
    av_coords = round(av_coords,0)
    if(my_img[matrix(av_coords,nrow=1)]!=unique_v[x]){
      dist = apply(t(t(my_coords)-av_coords)^2,1,sum)
      av_coords = my_coords[which.min(dist),]
    }
    out = av_coords + c(-92, -126, -72)
    return(out)
  }

  coords = lapply(1:length(unique_v),get.coords,unique_v,my_img)
  coords = do.call(rbind,coords)

  return(coords)
}
