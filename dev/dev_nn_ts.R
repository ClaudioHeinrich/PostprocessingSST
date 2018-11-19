rm(list = ls())

library(PostProcessing)
library(data.table)
library(keras)
library(tensorflow)

standard_mlp = function(X,
                        Y,
                        X_test,
                        hidden_config = 32,
                        epochs = 5,
                        batch_size = 1e3,
                        validate_prop = .1,
                        lambda = 0.01,
                        depth = 1,
                        verbose = 2){
    
    p_x = dim(X)[2]
    N = dim(X)[1]
    batch_size = 1e3
    
    ##if(is.null(batch_size)){batch_size = round(dim(X)[1]/1e2)}
    ##if(batch_size < 1)batch_size = 1
    inputs = layer_input(shape = p_x)
    predictions = inputs %>%
        layer_dense(units = hidden_config, activation = "relu",
                    kernel_regularizer = regularizer_l1(lambda),
                    use_bias = TRUE) %>%
        layer_dense(units = 1,use_bias = TRUE)
    
    model = keras_model(inputs = inputs, outputs = predictions)
    model %>% compile(loss = "MSE",
                      optimizer = "adam")
    verbose = 2
    history  = model %>% fit(as.matrix(X),Y, epochs = epochs, validation_split = .1,
                             batch_size = batch_size,
                             verbose = verbose)

    Y_test = model %>% predict(as.matrix(X_test))

    return(Y_test)
}

setwd("~/PostClimDataNoBackup/SFE/")
load("./FcNov2018/ts_hindcast_slimmed.RData")

for(m in 2:12){
    DT_final[,paste0("m",m) := (month == m) * 1]
}
nms_X = c("climatology", "ecmwf_anamoly", "norcpm_anamoly", paste0("m",2:12))
for(y in 2006:2017){
    print(y)
    DT_train = DT_final[year < y & year != min(year)]
    Y = DT_train[,obs_anamoly]
    X = DT_train[,.SD, .SDcols = nms_X]
    X_test = DT_final[year == y, .SD, .SDcols = nms_X]
    DT_final[year == y, p_hat := climatology + standard_mlp(X,Y,X_test, epochs = 10, lambda = .01, hidden_config = 8)]

    mod = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    DT_final[year == y, p_lm := climatology + predict(mod, newdata = DT_final[year == y])]
}

gg = function(a,b){return( sqrt(mean( (a - b)^2)))}

Lat_nordic = c(55,80)
Lon_nordic = c(5,30)
Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(p_hat),
                 lapply(.SD,gg,obs_erai_ts),
                 month,
                 .SDcols = c("climatology", "p_hat","p_lm")]

