# This is a script to generate the training dataset

# Varying parameters:
# sigma (volatilites)
# cor
# r
# K
# weights

# With default settings, each run of 2-stock FDS, returns 10,000 combinations

library(tidyverse)

source("fds_2_stock_basket.R")

# tst = fds_2s(K = 75, weights = c(.3, .7))

set.seed(420)

s1s = runif(2, .01, .99)

s2s = runif(2, .01, .99)

# https://www.sciencedirect.com/science/article/abs/pii/S0304405X10000437
cors = rnorm(3, .237, .093)

rs = runif(2, .01, .06)

Ks = runif(3, 10, 80)

weights = runif(3, .01, .99)

data = tibble(s1 = NA, s2 = NA, cor = NA, r = NA, K = NA, w1 = NA, w2 = NA, Stock_1_p = NA, Stock_2_p = NA, V = NA)

t1 = proc.time()

for (s1 in s1s) {
    for (s2 in s2s){
        for (c in cors) {
            for (r in rs) {
                for (strike in Ks) {
                    for (w1 in weights) {
                        ret = fds_2s(sigma = c(s1, s2), cor = c, r = r, K = strike, weights = c(w1, 1 - w1))
                        data = bind_rows(data, ret)
                    }
                }
            }
        }
    }
}

t2 = proc.time()


