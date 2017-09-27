##------- Cpp function equivalent ------
## First read the source code from a separate file
src <- paste(readLines("./example2.cpp"), collapse="\n")
## Now hand the source to cppFunction for compilation
cppFunction(src)
##-----------------------------------------
