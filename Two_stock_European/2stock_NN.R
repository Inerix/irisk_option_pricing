#####################
# @ Params
# sig: volatilty
# cor: Correlation between underlying
# r: Risk free interest rate
# K: Strike
# S: stock price
# w: weight
# weights: Vector representing weights of the stocks in the basket
#####################

library(devtools)
library(reticulate)
library(keras)
library(tensorflow)
library(tidyverse)

source("Two_stock_European/fds_2_stock_basket.R")

# data setup
data_2s = read.csv("data/two_stock.csv")
colnames(data_2s) = c("sig1", "sig2", "cor", "r", "K", "w1", "w2", "S1", "S2", "V")
n_2s = nrow(data_2s)
trn_2s_idx = sample(seq(from = 1, to = n_2s, by = 1), 0.8*n_2s)

# train-test split, removing w2 because of linear redundancy (w1 + w2 = 1)
data_2s_trn = as.matrix(data_2s[trn_2s_idx, -which(colnames(data_2s) == "w2")])
data_2s_tst = as.matrix(data_2s[-trn_2s_idx, -which(colnames(data_2s) == "w2")])

basket_mod = keras_model_sequential() %>% 
    layer_dense(units = 64, activation = "relu", input_shape = (ncol(data_2s_trn) - 1) ) %>% 
    layer_dense(units = 1) %>% 
    
    compile(
        optimizer = "rmsprop",
        loss = "mse",
        metrics = c("mae")  
    )

callbacks <- list(callback_early_stopping(patience = 12))
t1 = proc.time()
set.seed(42)
basket_fit = basket_mod %>%
    fit(
        x = data_2s_trn[, 1:8],
        y = data_2s_trn[, 9],
        epochs = 500, 
        batch = 1024,  
        validation_split = 0.25,
        callbacks = callbacks,
        verbose = TRUE
    )
t2 = proc.time()
t2-t1
# user: 810.169s

# training error
trn_2s_pred = basket_mod %>% predict(as.matrix(data_2s_trn[,1:8]))
trn_2s_act = data_2s_trn[,9]
plot(trn_2s_pred, trn_2s_act, main = "", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

trn_2s_rmse  = sqrt(mean((trn_2s_pred - trn_2s_act)^2))
trn_2s_rmse
# 0.5973254

# testing error
tst_2s_pred = basket_mod %>% predict(data_2s_tst[,1:8])
tst_2s_act = data_2s_tst[, 9]
plot(tst_2s_pred, tst_2s_act, main = "", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

tst_2s_rmse  = sqrt(mean((tst_2s_pred - tst_2s_act)^2))
tst_2s_rmse
#0.5979251

#####
# extrapolation

set.seed(400)
oob_sig1 = runif(1, .01, .99)
oob_sig2 = runif(1, .01, .99)
# https://www.sciencedirect.com/science/article/abs/pii/S0304405X10000437
oob_cor = rnorm(1, .237, .093)
oob_r = runif(1, .01, .06)
oob_K = runif(1, 10, 80)
oob_w = runif(1, .01, .99)
oob_data = tibble(sig1 = NA, sig2 = NA, cor = NA, r = NA, K = NA, w1 = NA, w2 = NA, S1 = NA, S2 = NA, V = NA)

ret = fds_2s(sigma = c(oob_sig1, oob_sig2), cor = oob_cor, r = oob_r, K = oob_K, weights = c(oob_w, 1 - oob_w))
oob_data = as.matrix(ret[, -which(colnames(data_2s) == "w2")])


oob_pred = basket_mod %>% predict(oob_data[,1:8])
oob_act = oob_data[, 9]
plot(tst_2s_pred, tst_2s_act, main = "oob actual vs pred", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

oob_rmse  = sqrt(mean((oob_act - oob_pred)^2))
oob_rmse
# 3.185204




