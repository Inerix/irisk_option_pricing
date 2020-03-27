# @ Params
# r - risk-free interest rate
# time - Vector representing time slices
# K - Strike price of the stock
# S - Vector of all possible stock prices
# sigma - Volatility 
# D - dividend rate
# @ Return - matrix of American Call option prices

######


######

# For convergence
# let a = 1/2 * sigma ** 2 * stock ** 2
# delta t < ds ** 2 / (2 * a)


#fds_1s_A(r, time, K, s, sigma, D)


fds_1s_A = function(r, time, K, S, sigma, D) {
    r = 0.06
    nstep_time = 1600
    time = seq(0, 1, length.out = nstep_time)
    sigma = .03
    D = .05
    K = 75
    S = seq(0, 150, length.out = 200)
    
    
    l_t = length(time)
    dt = time[2]
    ds = S[2]
    
    V = array(rep(0, (length(S)) * (l_t)), dim = c(l_t, length(S)))
    # [time, price]
    
    # First, consider final condition
    V[1, ] = pmax(S - K, 0)
    #mean(V[1, ] == 0)
    
    # Then, consider Boundary condition
    V[, 1] = pmax(S[1] - K * exp(-r * (T - time)), 0)
    #mean(V[,1 ] == 0)
    
    V[, length(S)] = pmax(S[length(S)] - K * exp(-r * time), 0)
    #mean(V[,401 ] == 0)
    
    V_amer = V
    # Finite difference scheme
    for (t in 1:(length(time) - 1)) {
        for (i in 2:(length(S) - 1)) {
            V[t + 1, i] =
                V[t, i - 1] * ((1 / 2) * ((sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2)) - ((r - D) * S[i] * dt) / (2 * ds)) + 
                V[t, i] * (1 - (sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2) - r * dt) +
                V[t, i + 1] * (dt / 2) * (((sigma ** 2 * S[i] ** 2) / (ds ** 2)) + ((r - D) * S[i]) / (ds))
            V_amer[t + 1, i] = max( round((S[i] - K), 4), round(V[t + 1, i], 4))
        }
    }
    mean(is.nan(V))
    return(V_amer)

}

    
