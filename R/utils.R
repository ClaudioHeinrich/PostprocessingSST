invert_ym = function(YM)
{
  ## TODO: verify
  return(c(YM %% 12, YM - YM %% 12 ))
}


restrict_coord = function(lons,lats)
{
  ss = parse(text = paste0("(Lon >",lons[1],") & (Lon <",lons[2],")",
                           "& (Lat >", lats[1],") & (Lat <",lats[2],")"))
  return(ss)
}

closest_point_helper = function(j, dt_1,dt_2, N = 1,
                                return_distance = FALSE,
                                distance_haversine = TRUE)
{
  if(j %% 1e2 == 0)print(paste("On",j))
  if(distance_haversine)
  {
    a = geosphere::distHaversine(as.vector(dt_1[j,.(Lon,Lat)]),as.matrix(dt_2[,.(Lon,Lat)]))
  }else{
    a = sqrt( ( dt_1[j,Lon] - dt_2[,Lon])^2 + (dt_1[j,Lat] - dt_2[,Lat])^2)
  }
  point_match = order(a)[1:N]
  if(return_distance)
  {
    point_match = c(point_match, a[point_match])
  }

  return(point_match)
}
