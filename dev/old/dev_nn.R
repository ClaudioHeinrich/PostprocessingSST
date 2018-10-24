rm(list = ls())

library(data.table)
library(keras)
library(tensorflow)

load("~/sfe_dt.RData")


standard_mlp = function(X,
                        Y,
                        f_loss,
                        p_model = 1,
                        X_test = NULL,
                        hidden_config = 32,
                        epochs = 5,
                        batch_size = NULL,
                        validate_prop = .1,
                        depth = 1,
                        verbose = 2){
    
    p_x = dim(X)[2]
    p_y = dim(Y)[2]

    N = dim(X)[1]
    
    if(is.null(batch_size)){batch_size = round(dim(X)[1]/1e2)}
    if(batch_size < 1)batch_size = 1
    
    inputs = layer_input(shape = p_x)
    if(depth ==0){
        predictions = inputs %>%
            layer_dense(units = p_model,use_bias = TRUE)
    }
    if(depth == 1){
        predictions = inputs %>%
            layer_dense(units = hidden_config, activation = "relu",
                        kernel_regularizer = regularizer_l1(0.01),
                        use_bias = TRUE) %>%
            layer_dense(units = p_model,use_bias = TRUE)
    }
    
    model = keras_model(inputs = inputs, outputs = predictions)
    model %>% compile(loss = f_loss,
                      optimizer = "adam")

    history  = model %>% fit(X,Y, epochs = epochs, validation_split = .1,
                             batch_size = batch_size,
                             verbose = verbose)

    l = list(model = model)
    if(!is.null(X_test)){
        theta_test = model %>% predict(X_test)
        l$theta_test = theta_test
    }

    return(l)
}


f_loss = function(Y_true, Yhat){

    mu = Yhat[,1]
    sigma2 = tf$exp(Yhat[,2])

    cond1 = tf$greater_equal(sigma2, 10000.0)
    sigma2 = tf$where(cond1,sigma2 - sigma2 + 10000, sigma2)

    cond2 = tf$less_equal(sigma2, .001)
    sigma2 = tf$where(cond2,sigma2 - sigma2 + .001, sigma2)

    y = Y_true[,1]

    
    l = -(y - mu)^2/ (2 * sigma2) - 1 /2 * tf$log(2 * pi * sigma2)

    return(tf$reduce_sum(-l, axis = 0L))
}

lloss = function(Y_true, Yhat){

    mu = Yhat[,1]
    sigma2 = exp(Yhat[,2])

    y = Y_true[,1]

    l = -(y - mu)^2/ (2 * sigma2) - 1 /2 * log(2 * pi * sigma2)

    return(sum(-l))
}


##Y_true = cbind(X %*% rep(1,dim(X)[2]),1)


##nms_features = c("Ens_bar","Ens_sd","Bias_Est","SD_hat")
nms_features = c("Ens_bar","Bias_Est")
nms_obs = "SST_bar"

S_curr = matrix(NA,12,2)
S_nn = matrix(NA,12,2)


R_all = NULL
for(m in 1:12){
    print(paste("M",m))
    YM_m = head(DT[year == 2010 &  month == m, YM],1)
    DT_m = DT[!is.na(Ens1) & YM <= YM_m]
    DT_m_train = DT_m[YM < YM_m & month == m]
    DT_m_test = DT_m[YM == YM_m]
    
    X_train = as.matrix(DT_m_train[, .SD, .SDcols = nms_features])
    Y_train = as.matrix(DT_m_train[, .SD, .SDcols = nms_obs])
    
    X_test = as.matrix(DT_m_test[, .SD, .SDcols = nms_features])
    
    
    l = list()
    X_train_sd = X_train
    Y_train_sd = Y_train
    mus = colMeans(X_train_sd)
    sds = apply(X_train_sd,2,"sd")
    
    for(j in 1:dim(X_train_sd)[2]){
        X_train_sd[,j] = (X_train_sd[,j] - mus[j]) / sds[j]
        X_test[,j] = (X_test[,j] - mus[j]) / sds[j]
    }
    
    Y_train_sd = (Y_train_sd - mean(Y_train_sd)) / sd(Y_train_sd)
    for(j in 1:5){
        print(paste("Attempt",j))
        l[[j]] = standard_mlp(X_train_sd,
                              Y_train_sd,
                              f_loss,
                              p_model = 2,
                              X_test = X_test,
                              hidden_config = 32,
                              batch_size = 1e3,
                              epochs = 20,
                              depth=1,
                              verbose = 2)
        
    }
    
    Y_test = as.matrix(DT_m_test[,.SD,.SDcols = nms_obs])
    S_nn = NULL
    Hat_nn = NULL
    Hat_simple = X_test[,1] * sds[1] + mus[1] + X_test[,2] * sds[2] + mus[2]
    for(j in 1:5){
        Hat_nn = cbind(Hat_nn, l[[j]]$theta_test[,1] * sd(Y_train) + mean(Y_train))
        S_nn[j] = sqrt(mean((Y_test - Hat_nn[,j])^2))
    }
    S_simple = sqrt(mean((Y_test - Hat_simple)^2))

    R_all = rbind(R_all, c(S_simple, S_nn))
}
    


