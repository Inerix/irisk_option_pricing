rm(list = ls())

source("Two_stock_European/fds_1_stock.R")

#####################
# @ Params
# sigma: Array of volatilties
# cor: Correlation between underlying
# r: Risk free interest rate
# K: Strike
# nstep_price: Number of price slices
# weights: Vector representing weights of the stocks in the basket
#####################

# one stock case
# a = 1/2 * \sigma ** 2 * S ** 2
# b = r * \sigma
# delta_s <= 2a/abs(b)

# fds_2s(sigma = c(s1, s2), cor, r, K, nstep_price, weights)

# initial values
fds_2s = function(sigma = c(.3, .3), cor = .3, r = .01, K, nstep_price = 100, weights = c(.5, .5)){
    s1 = sigma[1]
    s2 = sigma[2]
    w1 = weights[1]
    w2 = weights[2]
    
    # for 1 - stock case, convergence is nstep_time >= (sigma ** 2 * I ** 2)
    # Condition can be found on Wilmott pg. 1215
    T = 1
    nstep_time = (max(sigma) * nstep_price) ** 2 #900
    delta_t = T/nstep_time
    
    # Time vector increments
    time = seq(0, T, delta_t)
    
    V = array(rep(0, (nstep_price + 1) * (nstep_price + 1) * (nstep_time + 1)),
                       dim = c(nstep_time + 1, nstep_price + 1, nstep_price + 1))
    # dimensions are [time, s1, s2]
    
    # 3 is an arbitrary number determined by Wilmott pg. 1202
    max_b1 = w1 * K * 3
    b_1 = seq(0, max_b1, max_b1/nstep_price)

    max_b2 = w2 * K * 3
    b_2 = seq(0, max_b2, length.out = nstep_price + 1)
    
    ds1 = b_1[2]
    ds2 = b_2[2]
    
    for (i in 1:length(b_1)) {
        for (j in 1:length(b_2)) {
            V[1, i, j] = max( b_1[i] * w1 + b_2[j] * w2 - K, 0)
        }
    }
    
    # set 4 boundary conditions
    
    # sets min S1
    V[, 1, ] = fds_1s(r, time, K - w1 * b_1[1], w2 * b_2, s2) # both of these are 0
    
    # sets min S2
    V[, , 1] = fds_1s(r, time, K - w2 * b_2[1], w1 * b_1, s1) # 
    
    # sets max S1
    V[, nstep_price + 1, ] = fds_1s(r, time, K - w1 * b_1[nstep_price + 1], w2 * b_2, s2)
    
    # sets max S2
    V[, , nstep_price + 1] = fds_1s(r, time, K - w2 * b_2[nstep_price + 1] , w1 * b_1, s1)
    
    # iterate through at each slice
    
    for (t in 1:nstep_time) {
        for (i in 2:nstep_price) {
            # iterate through stock price 1s
            s1_t = b_1[i]
            s1_thing = (delta_t * (s1 * w1 * s1_t) ** 2) / ds1 ^ 2
            r_thing_1 = (r * delta_t * s1_t) / (2 * ds1) 
            for (j in 2:nstep_price) {
                s2_t = b_2[j]
                s2_thing = (delta_t * (s2 * w2 * s2_t) ** 2) / ds2 ^ 2
                r_thing_2 = (r * delta_t * s2_t) / (2 * ds2)
                huge_thing = (s1 * s2 * cor * w1 * w2 * s1_t * s2_t * delta_t) / 
                                (8 * ds1 * ds2)
        
                V[t + 1, i, j] = 
                    V[t, i, j] * (1 - s1_thing - s2_thing - r * delta_t) + 
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
    
    ret_dim_1 = floor(length(b_1) * w1 * 3 / 4)
    ret_dim_2 = floor(length(b_2) * w2  * 3 / 4)
    
    V = V[length(time), ret_dim_1:length(b_1), ret_dim_2:length(b_2)]
    
    ret_tibble = tibble(sig1 = rep(s1, length(V)), sig2 = rep(s2, length(V)), cor = rep(cor, length(V)), 
                        r = rep(r, length(V)), K = rep(K, length(V)), w1 = rep(w1, length(V)), w2 = rep(w2, length(V)),
                        S1 = rep(0, length(V)), S2 = rep(0, length(V)), 
                        V = rep(0, length(V)))
    alt_len = dim(V)[1]
    ret_tibble$V = c(V)
    ret_tibble$S1 = rep_len(b_1[ret_dim_1:length(b_1)], length.out = length(V))
    ret_tibble$S2 = unlist(lapply(b_2[ret_dim_2:length(b_2)], 
                                         function(price){
                                             rep(price, alt_len)
                                         }))
    
    return(ret_tibble)
}

