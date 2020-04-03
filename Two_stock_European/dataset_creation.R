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

set.seed(420)

sig1s = runif(2, .01, .99)        # 0.6034283, 0.9608715
sig2s = runif(2, .01, .99)        # 0.1809654, 0.4762730
# https://www.sciencedirect.com/science/article/abs/pii/S0304405X10000437
cors = rnorm(3, .237, .093)     # 0.2924472, 0.2079799, 0.2740121
rs = runif(2, .01, .06)         # 0.02516739, 0.04581248
Ks = runif(3, 10, 80)           # 30.33641, 23.73577, 32.58050
weights = runif(3, .01, .99)    # 0.8336193, 0.0115632, 0.1719843
data = tibble(sig1 = NA, sig2 = NA, cor = NA, r = NA, K = NA, w1 = NA, w2 = NA, S1 = NA, S2 = NA, V = NA)

t1 = proc.time()

for (s1 in sig1s) {
    for (s2 in sig2s){
        for (c in cors) {
            for (r in rs) {
                for (strike in Ks) {
                    for (w1 in weights) {
                        ret = fds_2s(sigma = c(sig1, sig2), cor = c, r = r, K = strike, weights = c(w1, 1 - w1))
                        data = bind_rows(data, ret)
                    }
                }
            }
        }
    }
}

t2 = proc.time()


