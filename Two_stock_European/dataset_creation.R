# This is a script to generate the training dataset

# Varying parameters:
# sigma (volatilites)
# cor
# r
# K
# weights

# With default settings, each run of 2-stock FDS, returns 10,000 combinations

library(tidyverse)

source("Two_stock_European/fds_2_stock_basket.R")

# tst = fds_2s(K = 75, weights = c(.3, .7))

dataset_size = 200000
n_combinations = 100
n_per_combo = dataset_size / n_combinations

set.seed(420)

sig1s = runif(n_combinations, .01, .99)
sig2s = runif(n_combinations, .01, .99)

# https://www.sciencedirect.com/science/article/abs/pii/S0304405X10000437
cors = rnorm(n_combinations, .237, .093)
rs = runif(n_combinations, .01, .06)
Ks = runif(n_combinations, 10, 80)
weights = runif(n_combinations, .05, .95)

data = tibble(sig1 = rep(NA, dataset_size), sig2 = rep(NA, dataset_size), cor = rep(NA, dataset_size), 
              r = rep(NA, dataset_size), K = rep(NA, dataset_size), w1 = rep(NA, dataset_size), 
              w2 = rep(NA, dataset_size), S1 = rep(NA, dataset_size), S2 = rep(NA, dataset_size), V = rep(NA, dataset_size))

t1 = proc.time()
for (ind in 0:(n_combinations - 1)) {
    ind_1_index = ind + 1
    data[(ind * n_per_combo + 1):((ind + 1) * n_per_combo ), ] = 
        fds_2s(sigma = c(sig1s[ind_1_index], sig2s[ind_1_index]), cor = cors[ind_1_index], r = rs[ind_1_index], 
               K = Ks[ind_1_index], weights = c(weights[ind_1_index], 1 - weights[ind_1_index]))
}

t2 = proc.time()


