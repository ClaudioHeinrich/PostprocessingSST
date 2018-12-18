rank1 = function(x){return(rank(x)[1])}
for(i in 1:length(q_print)){
    q = q_print[i]
    nm = paste0("p_below_",i)
    DT_forecast[,eval(nm) := apply(.SD,1,"rank1") / 1e2,
                .SDcols = c(paste0("q_",q), paste0("q_hat_",1:99))]    
}

save(DT_forecast,file = "./FcNov2018/Forecast_ts.RData")

