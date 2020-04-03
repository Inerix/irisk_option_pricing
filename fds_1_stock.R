# Adapting Biwen's code so that it is callable
# @ Params
# r - risk free interest rate
# time - Vector representing time slices
# strike - Strike price of the stock : K -  (w2 * s2) 
# S - Vector of all possible stock prices : (w1 * s1)
# sigma - Volatility (Black-Scholes)
#
# @ Return - 2-d matrix filled using FDS

######

# parameters taken from fds_2 issue plane
# strike = K - weights[1] * b_1[nstep_price + 1] # fix max stock 1
# b = weights[2] * b_2
# sigma = sigma[2]
######

# For convergence
# let a = 1/2 * sigma ** 2 * stock ** 2
# delta t < ds ** 2 / (2 * a)


#fds_1s(r, time, K - weights[1] * b_1[1], weights[2] * b_2, sigma[2])

fds_1s = function(r, time, strike, S, sigma) {
    l_t = length(time)
    dt = time[2]
    ds = S[2]
    
    V = array(rep(0, (length(S)) * (l_t)), dim = c(l_t, length(S)))
    # [time, price]
    
    # First, consider final condition
    V_c[1, ] = pmax(b - strike, 0)
    
    # Then, consider boundary condition
    V_c[, 1] = pmax(b[1] - strike * exp(-r * (T - time)), 0)
    
    V_c[, length(b)] = pmax(b[length(b)] - strike * exp(-r * time), 0)
    
    mat = array(rep(0, 3 * length(S)), c(3, length(S)))
    for (i in 1:length(S)) {
        mat[1, i] = (1 - (sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2) - r * dt)
        mat[2, i] = (dt / 2) * (((sigma ** 2 * S[i] ** 2) / (ds ** 2)) + (r * S[i]) / (ds))
        mat[3, i] = ((1 / 2) * ((sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2)) - (r * S[i] * dt) / (2 * ds))
    }
    # Finite difference scheme
    for (t in 1:(length(time) - 1)) {
        for (i in 2:(length(S) - 1)) {
            V[t + 1, i] =
                V[t, i] * (1 - (sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2) - r * dt) +
                V[t, i + 1] * (dt / 2) * (((sigma ** 2 * S[i] ** 2) / (ds ** 2)) + (r * S[i]) / (ds)) +
                V[t, i - 1] * ((1 / 2) * ((sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2)) - (r * S[i] * dt) / (2 * ds))
        }
    }
    return(V_c)
    
}
