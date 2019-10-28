rm(list = ls())

library(devtools,roxygen2)
library(data.table)

setwd('/nr/user/claudio/pkg/paper/PostprocessingSST')

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
