
##------ Plot ---------------
plot(dt_box[,.(Lat,Lon)], pch = ".", xlim = range(p_route[,1]),ylim = range(p_route[,2]))
map("world", add = TRUE)
points(nyc.loc[1],nyc.loc[2], pch = 19, col="blue")
points(bergen.loc[1], bergen.loc[2], pch = 19, col="blue")
points(p_route, pch = ".", col="blue")
segments(p_last[1],p_last[2], p_dest[1],p_dest[2], col="blue", lty =2)
p_last = p_dest


}


