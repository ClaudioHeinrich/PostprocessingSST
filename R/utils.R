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

closest_point_helper = function(j, dt_1,dt_2)
{
  if(j %% 1e2 == 0)print(paste("On",j))
  a = geosphere::distHaversine(as.vector(dt_1[j,.(Lon,Lat)]),as.matrix(dt_ens_grid[,.(Lon,Lat)]))
  point_match = order(a)[1]
  return(point_match)
}
