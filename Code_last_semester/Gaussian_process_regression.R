rm(list = ls())
finite_difference_scheme_upperbound = function(r, T, strike, S0, sigma, omega) {
    I = 200
    J = 400
    delta_t = 1 / 15000
    n = floor(T / delta_t)
    con = sqrt(2 * delta_t / (2 - r * delta_t))
    set.seed(321)
    delta_b = con / runif(1)
    time = seq(0, T, delta_t)
    b = seq(-I * delta_b, (J - I) * delta_b, delta_b)
    Num_B = length(b)
    V_c = matrix(NA, length(b), length(time))
    # First, consider final condition
    
    s1 = matrix(t(b) %o% sigma, length(b), length(sigma))  # make a matrix with in each row i: (b[i]*sigma1 b[i]*sigma2 b[i]*sigma3  b[i]*sigma4)
    s2 = matrix(t(rep(1, length(b))) %o% sigma, length(b), length(sigma))
    
    V_c[, 1] = pmax(rowSums(S0 * omega * exp((r - 0.5 * s2 ^ 2) * T + s1)) -
                        strike, 0)
    
    #for (i in 1:length(b)){
    #  V_c[i,1]=max(sum(S0*omega*exp((r-0.5*sigma^2)*T+sigma*b[i]))-strike,0)
    #}
    # Then, consider boundary condition
    
    V_c[1, ] = max(sum(S0 * omega * exp((r - 0.5 * sigma ^ 2) * T + sigma *
                                            b[1])) - strike, 0)
    V_c[J + 1, ] = max(sum(S0 * omega * exp((r - 0.5 * sigma ^ 2) * T + sigma *
                                                b[J + 1])) - strike, 0)
    
    # Finite difference scheme
    #for (k in 1:n){
    #  for (j in 2:J){
    #    V_c[j,k+1]=0.5*delta_t/delta_b^2*V_c[j-1,k]+(1-r*delta_t-delta_t/delta_b^2)*V_c[j,k]+0.5*delta_t/delta_b^2*V_c[j+1,k]  }
    #}
    
    V_New = c()
    V_Old = c()
    V_New = V_c[, 1]
    
    # Use the PDE in function of B.
    for (k in 1:n) {
        V_Old = V_New
        V_BB = (V_Old[3:(Num_B)] - 2 * V_Old[2:(Num_B - 1)] + V_Old[1:(Num_B -
                                                                           2)]) / delta_b ^ 2
        Theta = r * V_Old[2:(Num_B - 1)] - 0.5 * V_BB
        V_New[2:(Num_B - 1)] = V_Old[2:(Num_B - 1)] - delta_t * Theta
        V_New[1] = V_Old[1] * (1 - r * delta_t)
        #V_New[1]= max(sum(S0*omega*exp((r-0.5*sigma^2)*T+sigma*b[1]))-strike,0)
        V_New[J + 1] = 2 * V_New[Num_B - 1] - V_New[Num_B - 2]
        #V_New[J+1]=max(sum(S0*omega*exp((r-0.5*sigma^2)*T+sigma*b[J+1]))-strike,0)
    }
    
    # Get our aim upper bound
    #Aim_V_c=V_c[I+1,n+1]
    Aim_V_c = V_New[I + 1]
    return (Aim_V_c)
}

# set the parameters
r = 0.05
sigma = c(0.2, 0.3, 0.4, 0.2)
omega = c(0.25, 0.25, 0.25, 0.25)
S0 = c(100, 100, 100, 100)
S = sum(omega * S0)
K = 100

# determine the comonotonic basket option price using Monte Carlo Simulation to check if the finite difference gives a reasonable price
NSim = 10000
S_Sim = c()
Payoff = c()
Sims = rnorm(NSim)
for (i in 1:length(Sims)) {
    S_Sim[i] = sum(omega * S0 * exp((r - 0.5 * sigma ^ 2) * T + sigma * sqrt(T) *
                                        Sims[i]))
    Payoff[i] = max(S_Sim[i] - K, 0)
}
mean(S_Sim)
Call = exp(-r * T) * mean(Payoff)
Call

finite_difference_scheme_upperbound(r, T, K, S0, sigma, omega)


# First, compute new parameter nu_i (See Daniel's paper (Page 12, Equation 4.6))

finite_difference_scheme_lowerbound = function (r, T, strike, sigma, omega, rho) {
    sum = 0
    for (i in 1:4) {
        sum = sum + sum(sigma[i] * sigma)
    }
    nu = sum(S0 * rho * sigma) / sqrt((S[1] ^ 2 * rho * sum))
    # Get new sigma (See Daniel's paper(Page 5, Equation 2.7))
    new_sigma = nu * sigma
    y = finite_difference_scheme_upperbound(r, T, strike, S0, new_sigma, omega)
    return(y)
}

# Determine Upper and Lower bound
# set the parameters
r = 0.05
sigma = c(0.2, 0.3, 0.4, 0.2)
omega = c(0.25, 0.25, 0.25, 0.25)
S0 = c(100, 100, 100, 100)
S = sum(omega * S0)
rho = 0.3
K = seq(50, 150, 1)
Call_UB = c()
Call_LB = c()

for (i in 1:length(K)) {
    Call_UB[i] = finite_difference_scheme_upperbound(r, T, K[i], S0, new_sigma, omega)
    Call_LB[i] = finite_difference_scheme_lowerbound(r, T, K[i], S0, new_sigma, omega, rho)
}


# Produce sample path
trainsize = 100
trainr = runif(trainsize, min = 0.01, max = 0.04)
trainK = runif(trainsize, min = 80, max = 120)
trainT = runif(trainsize, min = 0.9, max = 1.1)

trainFDupper_bound = vector()
for (i in 1:trainsize) {
    r = trainr[i]
    K = trainK[i]
    Maturity = trainT[i]
    trainFDupper_bound[i] = finite_difference_scheme_upperbound(r, Maturity, K, S0, sigma, omega)
}

trainFDlower_bound = vector()
for (i in 1:trainsize) {
    r = trainr[i]
    K = trainK[i]
    Maturity = trainT[i]
    trainFDlower_bound[i] = finite_difference_scheme_lowerbound(r, Maturity, K)
}

library(GPFDA)
library(MASS)
trainX = cbind(trainr, trainT, trainK)
gprmd_UB = gpr(trainX, trainFDupper_bound)
gprmd_LB = gpr(trainX, trainFDlower_bound)

testsize = 100
set.seed(123)
testr = runif(testsize, min = 0.015, max = 0.035)
testK = runif(testsize, min = 90, max = 110)
testT = runif(testsize, min = 0.95, max = 1.05)
testX = cbind(testr, testT, testK)

t1 = Sys.time()
testFDupper_bound = vector()
for (i in 1:testsize) {
    r = testr[i]
    K = testK[i]
    Maturity = testT[i]
    testFDupper_bound[i] = finite_difference_scheme_upperbound(r, Maturity, K, S0, sigma, omega)
}
t2 = Sys.time()
FD_time = t2 - t1


testFDlower_bound = vector()
for (i in 1:testsize) {
    r = testr[i]
    K = testK[i]
    Maturity = testT[i]
    testFDlower_bound[i] = finite_difference_scheme_lowerbound(r, Maturity, K)
}

t3 = Sys.time()
test_pred_UB = gppredict(gprmd_UB, testX)[[1]]
t4 = Sys.time()
GPR_time = t4 - t3

test_pred_LB = gppredict(gprmd_LB, testX)[[1]]

x1 = test_pred_LB
y1 = testFDlower_bound
linearfit_LB = lm(y1 ~ x1)
plot(
    x1,
    y1,
    col = 'darkblue',
    main = 'Out-of-sample prediction for lower bound',
    xlab = 'GPR_LB',
    ylab = 'FD_LB'
)
abline(linearfit_LB$coefficients[1], linearfit_LB$coefficients[2], col =
           'red')

x2 = test_pred_UB
y2 = testFDupper_bound
linearfit_UB = lm(y2 ~ x2)
plot(
    x2,
    y2,
    col = 'darkorange',
    main = 'Out-of-sample prediction for upper bound',
    xlab = 'GPR_UB',
    ylab = 'FD_UB'
)
abline(linearfit_UB$coefficients[1], linearfit_UB$coefficients[2], col =
           'darkgreen')

AAE_LB = sum(abs(testFDlower_bound - test_pred_LB)) / testsize
MAE_LB = max(abs(testFDlower_bound - test_pred_LB))
plot(
    abs(testFDlower_bound - test_pred_LB),
    col = 'red',
    main = 'Error for lower bound',
    ylab = 'absolute error',
    xlab = 'the order of sample path'
)
AAE_UB = sum(abs(testFDupper_bound - test_pred_UB)) / testsize
MAE_UB = max(abs(testFDupper_bound - test_pred_UB))
plot(
    abs(testFDupper_bound - test_pred_UB),
    col = 'blue',
    main = 'Error for upper bound',
    ylab = 'absolute error',
    xlab = 'the order of sample path'
)


# Consider size influence
size_train = c(100, 500, 1000)
# generate train interest rate
number = max(size_train)
AllTrainr = runif(number, min = 0.01, max = 0.04)

# generate train strikes
AllTrainK = runif(number, min = 80, max = 120)

# generate train maturity
AllTrainT = runif(number, min = 0.9, max = 1.1)


AllFDtrain_UB = vector()
for (i in 1:number) {
    R = AllTrainr[i]
    MATURITY = AllTrainT[i]
    STRIKE = AllTrainK[i]
    AllFDtrain_UB[i] = finite_difference_scheme_upperbound(R, MATURITY, STRIKE)
}


AllFDtrain_LB = vector()
for (i in 1:number) {
    R = AllTrainr[i]
    MATURITY = AllTrainT[i]
    STRIKE = AllTrainK[i]
    AllFDtrain_LB[i] = finite_difference_scheme_lowerbound(R, MATURITY, STRIKE)
}
#Gaussian Process Regression
library(GPFDA)
DS_AAE_UB = vector()
for (i in 1:length(size_train)) {
    m = size_train[i]
    Data = cbind(AllTrainr[1:m], AllTrainT[1:m], AllTrainK[1:m])
    res_UB = AllFDtrain_UB[1:m]
    DS_AAE_UB[i] = sum(abs(gppredict(
        gpr(Data = Data, response = res_UB), testX
    )[[1]] - testFDupper_bound)) / testsize
}


DS_MAE_LB = vector()
for (i in 1:length(size_train)) {
    m = size_train[i]
    Data = cbind(AllTrainr[1:m], AllTrainT[1:m], AllTrainK[1:m])
    res_LB = AllFDtrain_LB[1:m]
    DS_AAE_LB[i] = sum(abs(gppredict(
        gpr(Data = Data, response = res_LB), testX
    )[[1]] - testFDlower_bound)) / testsize
}
