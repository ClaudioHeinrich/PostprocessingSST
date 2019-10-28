rm(list = ls())

library(devtools)
library(roxygen2)

#setwd('./')

ff = "./pp.sst"
if(!file.exists(ff)){
  dir.create('pp.sst')
}

# document

setwd('./pp.sst')
document()

setwd('..')
install('pp.sst')

library(pp.sst)
