rm(list = ls())

library(devtools,roxygen2)

setwd('/nr/user/claudio/pkg/paper_pkg/PostprocessingSST')

create('pp.sst')


# document

setwd('./pp.sst')
document()
