##########################################
# Scripts to run prior to this one
# | american-fds.R
##########################################
# @ Params
# n - number of combinations (# of observations)
# nn - n / 100 (100: number of S grid slices in FDS), all params except S has nn unique values
# r - risk-free interest rate
# time - Vector representing time slices
# K - Strike price
# S - Stock price (vector of size n)
# sigma - Volatility 
# D - dividend yield
# Convergence Condition: # time steps >= sigma^2 * i^2, i is number of stock steps
# @ Return - dataset of n observations

library(devtools)
library(reticulate)
library(keras)
library(tensorflow)
library(tidyverse)

source("American_options/american-fds.R")

obs_in_data = 60000
num_sets = obs_in_data/100

set.seed(42)
r = runif(num_sets, 0.01, 0.06)
sigma = runif(num_sets, 0.1, 0.6)
D = runif(num_sets, 0.01, 0.06)
K = runif(num_sets, 25, 75)

nstep_time = round(max(sigma), 1)^2 * 100^2   #3600

time = seq(0, 1, length.out = nstep_time)

data = data.frame("r" = rep(NA, obs_in_data), 
                  "sigma" = rep(NA, obs_in_data), 
                  "D" = rep(NA, obs_in_data), 
                  "K" = rep(NA, obs_in_data), 
                  "S" = rep(NA, obs_in_data),
                  "V" = rep(NA, obs_in_data))

t1 = proc.time()

for (k in 1:num_sets){
    # k goes from 1 to 100
    # i:j goes over 100 indices per each k
    # every set of 100 should have same r, sigma, D, K
    # within these sets of 100, we have a sequence of S and corresponding V
    i = 100*(k-1) + 1
    j = 100*k
    data$r[i:j] = rep(r[k], 100)
    data$sigma[i:j] = rep(sigma[k], 100)
    data$D[i:j] = rep(D[k], 100)
    data$K[i:j] = rep(K[k], 100)
    data$S[i:j] = seq(0, 3*K[k], length.out = 100)
    data$V[i:j] = fds_1s_A(time, r[k], K[k], sigma[k], D[k])
}

t2 = proc.time()
t2-t1
# user time: 541.604

set.seed(42)
trn_idx = sample(1:obs_in_data, 0.8 * obs_in_data)
data_trn = as.matrix(data[trn_idx,])
data_tst = as.matrix(data[-trn_idx,])

set.seed(42)
amer_mod = keras_model_sequential() %>% 
    layer_dense(units = 64, activation = "relu", input_shape = 5) %>% 
    layer_dense(units = 1, activation = "relu") %>% 
    compile(
        optimizer = "rmsprop",
        loss = "mse",
        metrics = c("mae")  
    )

t3 = proc.time()

callbacks <- list(callback_early_stopping(patience = 20))

set.seed(42)
amer_fit = amer_mod %>%
    fit(
        x = data_trn[, 1:5],
        y = data_trn[, 6],
        epochs = 100, 
        batch = 1024,  
        validation_split = 0.25,
        callbacks = callbacks,
        verbose = TRUE
    )
t4 = proc.time()
t4-t3

# training error
trn_pred = amer_mod %>% predict(data_trn[,1:5])
plot(trn_pred, data_trn[,6], main = "", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

trn_mae  = mean(abs((trn_pred - data_trn[,6])))
trn_mae
# 1.056866

# tsting error
tst_pred = amer_mod %>% predict(data_tst[,1:5])
plot(tst_pred, data_tst[,6], main = "test", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")d

tst_mae  = mean(abs((tst_pred - data_tst[,6])))
tst_mae
# 1.06437

write.csv(data, "American_options/data/american_data.csv")
