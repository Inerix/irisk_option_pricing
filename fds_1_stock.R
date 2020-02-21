# Adapting Biwen's code so that it is callable
# @ Params
# r - risk free interest rate
# time - Vector representing time slices
# strike - Strike price of the stock (w2 * s2) - K
# b - Vector of all possible stock prices (w1 * s1)
# sigma - Volatility (Black-Scholes)
#
# @ Return - 2-d matrix filled using FDS

fds_1s = function(r, time, strike, b, sigma, dt, db) {
    l_t = length(time)
    
    V_c = array(rep(0, (length(b)) * (l_t)), dim = c(l_t, length(b)))
    # [time, price]

    # First, consider final condition
    V_c[1, ] = max(b - strike, 0)
    
    # Then, consider boundary condition
    V_c[, 1] = max(b[1] - strike, 0)
    V_c[, length(b)] = max(b[length(b)] - strike, 0)
    
    # Finite difference scheme
    for (t in 1:(length(time) - 1)) {
        for (i in 2:(length(b) - 1)) {
            V_c[t + 1, i] = 
                V_c[t, i] * (1 - (sigma ^ 2 * b[i] ^ 2 * dt) / (db ^ 2) - r * dt) + 
                V_c[t, i + 1] * ((1 / 2) * (sigma ^ 2 * b[i] ^ 2 * dt) / (db ^ 2) + (r * b[i] * dt) / (2 * db)) + 
                V_c[t, i - 1] * ((1 / 2) * (sigma ^ 2 * b[i] ^ 2 * dt) / (db ^ 2) + (r * b[i] * dt) / (2 * db))
        }
    }

    return(V_c)
}