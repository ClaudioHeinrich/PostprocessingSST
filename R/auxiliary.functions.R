


#' truncation function for freezing level of sea water
#' 
#' @author Claudio Heinrich
#' 
#' @example trc(c(-2,3))
#' 
#' @export

trc = function (x){ 
  truncation.value = -1.769995
  x = truncation.value * (x < truncation.value) + x * (x >= truncation.value)
  return(x)}