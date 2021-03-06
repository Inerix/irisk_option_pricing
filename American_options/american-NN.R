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
    # i:j goes over 100 indices per each k
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

trn_idx = sample(1:obs_in_data, 0.8 * obs_in_data)
data_est = as.matrix(data[trn_idx,])
data_val = as.matrix(data[-trn_idx,])

amer_mod = keras_model_sequential() %>% 
    layer_dense(units = 64, activation = "relu", input_shape = 5) %>% 
    layer_dense(units = 1) %>% 
    compile(
        optimizer = "rmsprop",
        loss = "mse",
        metrics = c("mae")  
    )

t3 = proc.time()

callbacks <- list(callback_early_stopping(patience = 12))

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
est_pred = amer_mod %>% predict(data_est[,1:5])
plot(est_pred, data_est[,6], main = "", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

est_rmse  = sqrt(mean((est_pred - data_est[,6])^2))
est_rmse
# 0.7316702

# validation error
val_pred = amer_mod %>% predict(data_val[,1:5])
plot(val_pred, data_val[,6], main = "", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

val_rmse  = sqrt(mean((val_pred - data_val[,6])^2))
val_rmse
#0.736379


### extrapolation test sets

n_tst1 = 0.2 * obs_in_data
nn_tst1 = n_tst1 / 100

set.seed(42)
r_tst1 = runif(n_tst1, 0.035, 0.085)
sigma_tst1 = runif(n_tst1, 0.35, 0.85)
D_tst1 = runif(n_tst1, 0.035, 0.085)
K_tst1 = runif(n_tst1, 50, 100)

nstep_time1 = round(max(sigma_tst1), 1)^2 * 100^2   #6400
time1 = seq(0, 1, length.out = nstep_time1)

S_tst1 = c(rep(NA,n_tst1))
data_tst1 = data.frame("r" = rep(NA, n_tst1), 
                  "sigma" = rep(NA, n_tst1), 
                  "D" = rep(NA, n_tst1), 
                  "K" = rep(NA, n_tst1), 
                  "S" = rep(NA, n_tst1),
                  "V" = rep(NA, n_tst1))

for (k in 1:nn_tst1){
    # i:j goes over 100 indices per each k
    i = 100*(k-1) + 1
    j = 100*k
    print(k)
    print(i)
    print(j)
    data_tst1$r[i:j] = rep(r_tst1[k], 100)
    data_tst1$sigma[i:j] = rep(sigma_tst1[k], 100)
    data_tst1$D[i:j] = rep(D_tst1[k], 100)
    data_tst1$K[i:j] = rep(K_tst1[k], 100)
    data_tst1$S[i:j] = seq(0, 3*K_tst1[k], length.out = 100)
    data_tst1$V[i:j] = fds_1s_A(time1, r_tst1[k], K_tst1[k], sigma_tst1[k], D_tst1[k])
}

# mean(is.nan(data_tst1$V))

tst1_pred = amer_mod %>% predict(as.matrix(data_tst1[,1:5]))
tst1_act = as.matrix(data_tst1[,6])
plot(tst1_pred, tst1_act, main = "", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

tst1_rmse  = sqrt(mean((tst1_pred - tst1_act)^2))
tst1_rmse
#0.736379