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
