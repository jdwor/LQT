#' @keywords internal
#' @importFrom neurobase readnii writenii

util_get_coords=function(cfg){
  # load parcels
  my_img = readnii(cfg$parcel_path)

  # get unique values
  unique_v = unique(my_img);
  unique_v = unique_v[unique_v > 0]

  # loop through unique values and get coordinates in MNI
  coords = matrix(nrow=length(unique_v),ncol=3)
  for(i in 1:length(unique_v)){
    my_coords = which(my_img==unique_v[i],arr.ind=T)
    av_coords = apply(my_coords,2,mean)
    av_coords = round(av_coords,0)
    if(my_img[my_coords]!=unique_v[i]){
      dist = apply(t(t(my_coords)-av_coords)^2,1,sum)
      av_coords = my_coords[which.min(dist),]
    }
    coords[i,] = av_coords + c(-92, -126, -72);
  }

  return(coords)
}
