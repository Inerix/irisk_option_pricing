# Biwen's code

# we consider comonotonic situation.

# The upper bound

# Initial value
S0 = c(100, 100, 100, 100)

# Maturity
T = 1

# Choose some integer I and J.
I = 200
J = 400

# Risk-free interest rate
r = 0.01

# Consider each time step
delta_t = 1 / 15000

# The number of steps
n = floor(T / delta_t)

# Sigma for black-scholes model (Volatility)
sigma = c(0.2, 0.3, 0.4, 0.2)

# The weight
omega = c(0.25, 0.25, 0.25, 0.25)

# The weighted stock price
S = omega * S0

# Strike price
strike = 100

# correlation
cor = 0.3

# Consider convergence condition.
con = sqrt(2 * delta_t / (2 - r * delta_t))

set.seed(321)
delta_b = con / runif(1)

# The time grids and shifted brownian motion grids
time = seq(0, T, delta_t)

# b is underlying price
b = seq(-I * delta_b, (J - I) * delta_b, delta_b)

## Start to generate upper bound
V_c = matrix(NA, length(b), length(time))

# First, consider final condition
for (i in 1:length(b)) {
    V_c[i, 1] = max(sum(S0 * omega * exp((r - 0.5 * sigma ^ 2) * T + sigma * b[i])) - strike, 0)
}

# Then, consider boundary condition
V_c[1,] = max(sum(S0 * omega * exp((r - 0.5 * sigma ^ 2) * T + sigma * b[1])) -
                  strike, 0)
V_c[J + 1,] = max(sum(S0 * omega * exp((r - 0.5 * sigma ^ 2) * T + sigma * b[J +
                                                                                 1])) - strike, 0)

# Finite difference scheme
for (k in 1:n) {
    for (j in 2:J) {
        V_c[j, k + 1] = 0.5 * delta_t / delta_b ^ 2 * V_c[j - 1, k] + (1 - r * delta_t -
                                                                           delta_t / delta_b ^ 2) * V_c[j, k] + 0.5 * delta_t / delta_b ^ 2 * V_c[j +
                                                                                                                                                      1, k]
    }
}

# Get our aim
Aim_V_c = V_c[I + 1, n + 1]

## Start to generate lower bound
# First, compute new parameter nu_i (See Daniel's paper (Page 12, Equation 4.6))
sum = 0
for (i in 1:4) {
    sum = sum + sum(sigma[i] * sigma)
}
nu = sum(S * cor * sigma) / sqrt((S[1] ^ 2 * cor * sum))
# Get new sigma (See Daniel's paper(Page 5, Equation 2.7))
new_sigma = nu * sigma
# Start to generate lower bound
V_l = matrix(NA, length(b), length(time))
# Consider final condition
for (i in 1:length(b)) {
    V_l[i, 1] = max(sum(S0 * omega * exp((
        r - 0.5 * new_sigma ^ 2
    ) * T + new_sigma * b[i])) - strike, 0)
}
# Consider boundary conditions
V_l[1, ] = max(sum(S0 * omega * exp((r - 0.5 * new_sigma ^ 2) * T + new_sigma *
                                        b[1])) - strike, 0)
V_l[J + 1, ] = max(sum(S0 * omega * exp((r - 0.5 * new_sigma ^ 2) * T +
                                            new_sigma * b[J + 1])) - strike, 0)

# Finite difference schemes
for (k in 1:n) {
    for (j in 2:J) {
        V_l[j, k + 1] = 0.5 * delta_t / delta_b ^ 2 * V_l[j - 1, k] + (1 - r * delta_t -
                                                                           delta_t / delta_b ^ 2) * V_l[j, k] + 0.5 * delta_t / delta_b ^ 2 * V_l[j +
                                                                                                                                                      1, k]
    }
}

# Get lower bound
Aim_V_l = V_l[I + 1, n + 1]
