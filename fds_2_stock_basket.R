rm(list = ls())
library(doParallel)
library(foreach)
source("fds_1_stock.R")

# in previous FDS, considered stock price of 1 commonotonic price over 400 steps and
# and considered price over 15,000 steps
# in our case this lead to 2,400,000,000 obsercations

# starting off with 100 * 100 * 3750 (37,500,000)

#####################
# @ Params
# nstock: Number of stock price slices
# weights: vector of weights (s1, s2)
# S_1: Initial price S1 (May or may not be necessary)
# S_2: Initial price S2 (May or may not be necessary)
#####################

# initial values

T = 1
K = 75
S_1 = 100
S_2 = 50
sigma = c(.2, .3) 
cor = 0.3
weights = c(.5, .5)

# Risk-free interest rate
r = 0.01

nstep_price = 100 # delta_s <= sigma ** 2 * S / r

# for 1 - stock case, convergence is delta_t <= (sigma ** 2 * I ** 2)
# Condition can be found on Wilmott pg. 1258
#delta_t = 1/((min(sigma) ** 2) * (nstep_price ** 2))
#nstep_time = 1/delta_t # 1600
nstep_time = 64000
delta_t = 1/nstep_time

t1 = proc.time()

V = array(rep(0, (nstep_price + 1) * (nstep_price + 1) * (nstep_time + 1)),
                   dim = c(nstep_time + 1, nstep_price + 1, nstep_price + 1))
# dimensions are [time, s1, s2]

# delta b_1 needs to be set so that the price range goes from 0 - 2 * stock price

# 3 is an arbitrary number determined by Wilmott pg. 1202
# once this is converted into a function this will me 
#       max(weights[1] * K * 3, S_1)
max_b1 = weights[1] * K * 3

b_1 = seq(0, max_b1, max_b1/nstep_price)

delta_S1 = b_1[2]

# once this is converted into a function this will me 
#       max(weights[2] * K * 3, S_2)
max_b2 = weights[2] * K * 3

b_2 = seq(0, max_b2, length.out = nstep_price + 1)

delta_S2 = b_2[2]

# For time to converge
# Wilmott pg. 1253
# a = coefficient on d^2v/dS1^2
a = .5 * (sigma[1] ** 2) * (weights[1] ** 2) * (S_1 ** 2)

# d = coefficient on d^2v/dS2^2

d = .5 * (sigma[2] ** 2) * (weights[2] ** 2) * (S_2 ** 2)

# Time vector increments
time = seq(0, T, delta_t)

# foreach (i = 1:length(b_1)) %dopar%{
for (i in 1:length(b_1)) {
    for (j in 1:length(b_2)) {
        V[1, i, j] = max( b_1[i] * weights[1] + b_2[j] * weights[2] - K, 0)
    }
}

# set 4 boundary conditions

# sets min S1
V[, 1, ] = fds_1s(r, time, K - weights[1] * b_1[1], weights[2] * b_2, sigma[2]) # both of these are 0
#mean(is.nan(V[, 1, ]))

# sets min S2
V[, , 1] = fds_1s(r, time, K - weights[2] * b_2[1], weights[1] * b_1, sigma[1]) # 
#mean(is.nan(V[, , 1]))

# sets max S1
V[, nstep_price + 1, ] = fds_1s(r, time, K - weights[1] * b_1[nstep_price + 1], weights[2] * b_2, sigma[2])
#mean(is.nan(V[, nstep_price + 1, ]))

# sets max S2
V[, , nstep_price + 1] = fds_1s(r, time, K - weights[2] * b_2[nstep_price + 1] , weights[1] * b_1, sigma[1])
#mean(is.nan(V[, , nstep_price + 1 ]))

# iterate through at each slice

ds1 = b_1[2]
ds2 = b_2[2]

#cl = makeCluster(2, type = "FORK")
#registerDoParallel(cores = 3)

#getDoParWorkers()

for (t in 1:nstep_time) {
    for (i in 2:nstep_price) {
        # iterate through stock price 1s
        s1_t = b_1[i]
        s1_thing = (delta_t * (sigma[1] * weights[1] * s1_t) ** 2) / ds1 ^ 2
        r_thing_1 = (r * delta_t * s1_t) / (2 * ds1) 
        for (j in 2:nstep_price) {
        #V[t + 1, i, 2:nstep_price] = foreach (j = 2:nstep_price, .combine = c) %dopar%{
            s2_t = b_2[j]
            s2_thing = (delta_t * (sigma[2] * weights[2] * s2_t) ** 2) / ds2 ^ 2
            r_thing_2 = (r * delta_t * s2_t) / (2 * ds2)
            huge_thing = (sigma[1] * sigma[2] * cor * weights[1] * weights[2] * s1_t * s2_t * delta_t) / 
                            (8 * ds1 * ds2)
            #return(
            V[t + 1, i, j] = 
                V[t, i, j] * (1 - s1_thing - s2_thing - r * delta_t) + 
                V[t, i + 1, j] * (r_thing_1 + s1_thing / 2) + 
                V[t, i - 1, j] * (s1_thing / 2 - r_thing_1) + 
                V[t, i, j + 1] * (r_thing_2 + s2_thing/2) + 
                V[t, i, j - 1] * (s2_thing / 2 - r_thing_2) + 
                V[t, i + 1, j + 1] * huge_thing - 
                V[t, i - 1, j + 1] * huge_thing - 
                V[t, i + 1, j - 1] * huge_thing +
                V[t, i - 1, j - 1] * huge_thing#)
        }    
    }
}

#stopImplicitCluster()

# total nan 

mean(is.nan(V))

t2 = proc.time()

#write.csv2(V, file = "fds_out.csv")

