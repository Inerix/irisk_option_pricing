rm(list = ls())
library(doMC)
library(foreach)
source("fds_1_stock.R")

# in previous FDS, considered stock price of 1 commonotonic price over 400 steps and
# and considered price over 15,000 steps
# in our case this lead to 2,400,000,000 obsercations

# starting off with 100 * 100 * 3750 (37,500,000)

nstep_time = 3750
nstep_price = 100

V = array(rep(0, (nstep_price + 1) * (nstep_price + 1) * (nstep_time + 1)), 
                     dim = c(nstep_time + 1, nstep_price + 1, nstep_price + 1))
# dimensions are [time, s1, s2]

# initial values

T = 1
delta_t = T/nstep_time
K = 75
S_1 = 100
S_2 = 50
sigma = c(.2, .3) 
cor = 0.3
weights = c(.5, .5)

# Risk-free interest rate
r = 0.01

# Time vector increments
time = seq(0, T, delta_t)

# delta_b_1 represents the total+- in stock 1
delta_b_1 = 30

b_1 = seq(S_1 - delta_b_1/2, S_1 + delta_b_1/2, length.out = nstep_price + 1)

# delta_b_2 represents the +- in stock 2
delta_b_2 = 10

b_2 = seq(S_2 - delta_b_2/2, S_2 + delta_b_2/2, length.out = nstep_price + 1)

# insert final condition

# cl <- makeCluster(2)
# registerDoParallel(cl)

# foreach (i = 1:length(b_1)) %dopar%{
for (i in 1:length(b_1)) {
    for (j in 1:length(b_2)) {
        V[1, i, j] = max( b_1[i] * weights[1] + b_2[j] * weights[2] - K, 0)
    }
}

# set 4 boundary conditions

# sets min S1
V[, 1, ] = fds_1s(r, time, K - weights[1] * b_1[1], weights[2] * b_2, sigma[2])

mean(is.nan(V[, 1, ]))

# sets min S2
V[, , 1] = fds_1s(r, time, K - weights[2] * b_2[1], weights[1] * b_1, sigma[1])

mean(is.nan(V[, , 1]))

# sets max S1
V[, nstep_price + 1, ] = fds_1s(r, time, K - weights[1] * b_1[nstep_price + 1], weights[2] * b_2, sigma[2])

mean(is.nan(V[, nstep_price + 1, ]))

# sets max S2
V[, , nstep_price + 1] = fds_1s(r, time, K - weights[2] * b_2[nstep_price + 1] , weights[1] * b_1, sigma[1])

mean(is.nan(V[, , nstep_price + 1]))

# iterate through at each slice

ds1 = delta_b_1 
ds2 = delta_b_2

for (t in 1:nstep_time){
    for (i in 2:nstep_price) {
        # iterate through stock price 1s
        s1_t = b_1[i]
        s1_thing = ((sigma[1] * weights[1] * s1_t) ^ 2) / ds1 ^ 2
        r_thing_1 = (r * delta_t * s1_t) / (2 * ds1) 
        for (j in 2:nstep_price) {
        #foreach (i = 2:nstep_price) %dopar%{
            # can do this step in parllel
            s2_t = b_2[j]
            s2_thing = ((sigma[2] * weights[2] * s2_t) ^ 2) / ds2 ^ 2
            r_thing_2 = (r * delta_t * s2_t) / (2 * ds2)
            huge_thing = (sigma[1] * sigma[2] * cor * weights[1] * weights[2] * s1_t * s2_t * delta_t) / 
                            (8 * delta_b_1 * delta_b_2)
            V[t + 1, i, j] = 
                V[t, i, j] * (1 - s1_thing + s2_thing - r * delta_t) + 
                V[t, i + 1, j] * (r_thing_1 + s1_thing / 2) + 
                V[t, i - 1, j] * (s1_thing / 2 - r_thing_1) + 
                V[t, i, j + 1] * (r_thing_2 + s2_thing/2) + 
                V[t, i, j - 1] * (s2_thing / 2 - r_thing_2) + 
                V[t, i + 1, j + 1] * huge_thing - 
                V[t, i - 1, j + 1] * huge_thing - 
                V[t, i + 1, j - 1] * huge_thing +
                V[t, i - 1, j - 1] * huge_thing
        }    
    }
}

# total nan 

mean(is.nan(V))

