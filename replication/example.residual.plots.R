

#---- get example residual plots for the years and months specified in Yvec and Mvec ----

dt = load_combined_wide(bias = TRUE)
dt[,"res":= SST_hat - SST_bar]

Yvec = c(2005,2007,2009)
Mvec = 9

#--- get maximum range for uniform plotting scale ----

rr = c()
for(Y in Yvec){
  for(M in Mvec){
    rr = c(rr,range(dt[year == Y & month == M,res],na.rm = TRUE))
  }
}
  

rr = range(rr)


for(Y in Yvec){
  for(M in Mvec){
    
    plot_system(Y = Y, M=M, type = "res", print_figs = FALSE, plot_title = paste0("Residual 0",M," / ",Y),rr = rr )
  }
}
