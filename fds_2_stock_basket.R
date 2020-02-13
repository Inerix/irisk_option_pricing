library(doMC)
library(foreach)

rm(list = ls())
# in previous FDS, considered stock price of 1 commonotonic price over 400 steps and
# and considered price over 15,000 steps
# in our case this lead to 2,400,000,000 obsercations

# starting off with 100 * 100 * 3750 (37500000)

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

# convergence
con = sqrt(2 * delta_t / (2 - r * delta_t))

set.seed(321)
delta_b_1 = con / runif(1)

time = seq(0, T, delta_t)

b_1 = seq(-nstep_price/2 * delta_b_1, nstep_price/2 * delta_b_1, delta_b_1)

delta_b_2 = con / runif(1)

b_2 = seq(-nstep_price/2 * delta_b_2, nstep_price/2 * delta_b_2, delta_b_2)

# insert final condition

# cl <- makeCluster(2)
# registerDoParallel(cl)

# foreach (i = 1:length(b_1)) %dopar%{
for (i in 1:length(b_1)) {
    for (j in 1:length(b_2)) {
        V[1, i, j] = max( S_1 * weights[1] * exp((r - 0.5 * sigma[1] ^ 2) * T + sigma[1] * b_1[i]) +
                          S_2 * weights[2] * exp((r - 0.5 * sigma[2] ^ 2) * T + sigma[2] * b_2[j]) - K, 0)
    }
}

# set 4 boundary conditions

# sets min S1
for (j in 1:length(b_2)) {
    val_j = max( S_1 * weights[1] * exp((r - 0.5 * sigma[1] ^ 2) * T + sigma[1] * b_1[1]) +
                     S_2 * weights[2] * exp((r - 0.5 * sigma[2] ^ 2) * T + sigma[2] * b_2[j]) - K, 0)
    for (t in 1:length(time)) {
        V[t, 1, j] = val_j
    }
}

# sets min S2
for (i in 1:length(b_1)) {
    val_i = max( S_1 * weights[1] * exp((r - 0.5 * sigma[1] ^ 2) * T + sigma[1] * b_1[i]) +
                     S_2 * weights[2] * exp((r - 0.5 * sigma[2] ^ 2) * T + sigma[2] * b_2[1]) - K, 0)
    for (t in 1:length(time)) {
        V[t, i, 1] = val_i
    }
}

# sets max S1
for (j in 1:length(b_2)) {
    val_j = max( S_1 * weights[1] * exp((r - 0.5 * sigma[1] ^ 2) * T + sigma[1] * b_1[nstep_price + 1]) +
                     S_2 * weights[2] * exp((r - 0.5 * sigma[2] ^ 2) * T + sigma[2] * b_2[j]) - K, 0)
    for (t in 1:length(time)) {
        V[t, nstep_price + 1, j] = val_j
    }
}

# sets max S2
for (i in 1:length(b_1)) {
    val_i = max( S_1 * weights[1] * exp((r - 0.5 * sigma[1] ^ 2) * T + sigma[1] * b_1[i]) +
                     S_2 * weights[2] * exp((r - 0.5 * sigma[2] ^ 2) * T + sigma[2] * b_2[nstep_price + 1]) - K, 0)
    for (t in 1:length(time)) {
        V[t, i, nstep_price + 1] = val_i
    }
}

# iterate through at each slice

ds1 = delta_b_1 
ds2 = delta_b_2

for (t in 2:length(time)){
    for (i in 1:length(b_1)) {
        s1_t = S_1 * exp((r - 0.5 * sigma[1] ^ 2) * (T - time[t]) + sigma[1] * b_1[i])

        s1_thing = ((sigma[1] * weights[1] * s1_t) ^ 2) / 
        for (j in 1:length(b_2)) {
            # can do this step in parllel
            s2_t = S_2 * exp((r - 0.5 * sigma[2] ^ 2) * (T - time[t]) + sigma[2] * b_2[j])
            V[t, i, j] = 
                V[t-1, i, j] * (1 - )
        }    
    }
}

