invert_ym = function(YM)
{
  ## TODO: verify
  return(c(YM %% 12, YM - YM %% 12 ))
}
