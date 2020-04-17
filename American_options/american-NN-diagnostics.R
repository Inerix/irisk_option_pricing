##########################################
# Scripts to run prior to this one
# | american-fds.R
# | american-NN.R

# remove data to clear up some memory
rm(data)
rm(data_trn)
rm(data_tst)
library(ggplot2)

################################################
# Diagnostics of American Neural Network Model
# || I: Explores how NN performs with each parameter
# || II: Explores extrapolation

################################################################################################
# I: Exploring Parameters
# 1) hold all but 1 of the 5 parameters constant
# 2) feed into our NN model to see how it performs
# 3) repeat for each parameter
################################################################################################

n_obs_ep = 10000
n_sets_ep = n_obs_ep/100

# non-constant parameters
set.seed(490)
r_ep = seq(from = 0.01, to = 0.06, length.out = n_sets_ep)
sigma_ep = seq(from = 0.1, to = 0.6, length.out = n_sets_ep)
D_ep = seq(from = 0.01, to = 0.06, length.out = n_sets_ep)
K_ep = seq(from = 25, 75, length.out = n_sets_ep)
data_ep_nc = data.frame("r" = r_ep, "sigma" = sigma_ep,"D" = D_ep,"K" = K_ep)

# constant parameters (median, index = 250)
cp_idx = floor(n_sets_ep/2)
cp = data.frame("r" = r_ep[cp_idx],
                "sigma" = sigma_ep[cp_idx],
                "D" = D_ep[cp_idx],
                "K" = K_ep[cp_idx])

nstep_time_ep = round(max(sigma_ep), 1)^2 * 100^2   #3600
time = seq(0, 1, length.out = nstep_time_ep)

####################################
# non-constant: r, risk-free rate
####################################

data_nc_r = data.frame("r" = r_ep, 
                       "sigma" = rep(NA, n_obs_ep),
                       "D" = rep(NA, n_obs_ep),
                       "K" = rep(NA, n_obs_ep),
                       "S" = rep(NA, n_obs_ep),
                       "V" = rep(NA, n_obs_ep))

t1 = proc.time()
for (k in 1:n_sets_ep){
    # k goes from 1 to 100
    # i:j goes over 100 indices per each k
    # every set of 100 should have same r, sigma, D, K
    # within these sets of 100, we have a sequence of S and corresponding V
    i = 100*(k-1) + 1
    j = 100*k
    data_nc_r$r[i:j] = rep(r_ep[k], 100)
    data_nc_r$sigma[i:j] = rep(cp$sigma, 100)
    data_nc_r$D[i:j] = rep(cp$D, 100)
    data_nc_r$K[i:j] = rep(cp$K, 100)
    data_nc_r$S[i:j] = seq(0, 3*cp$K, length.out = 100)
    data_nc_r$V[i:j] = fds_1s_A(time, r_ep[k], cp$K, cp$sigma, cp$D)
}
t2 = proc.time()
t2-t1

data_nc_r = as.matrix(data_nc_r)
r_pred = amer_mod %>% predict(data_nc_r[,1:5])
plot(r_pred, data_nc_r[,6], main = "r, non-constant", xlab = "actual (V)", ylab = "predicted (V)", cex = 0.5)
abline(0, 1, col = "red")

r_mae  = mean(abs((r_pred - data_nc_r[,6])))
r_mae
# 0.957807

# get median S for each set of 100
nc_r = data.frame("r" = rep(NA, n_sets_ep),
                  "sigma" = rep(NA, n_sets_ep),
                  "D" = rep(NA, n_sets_ep),
                  "K" = rep(NA, n_sets_ep),
                  "S" = rep(NA, n_sets_ep),
                  "V" = rep(NA, n_sets_ep),
                  "V_pred" = rep(NA, n_sets_ep))

for (k in 1:n_sets_ep){
    # k goes from 1 to 100
    # i = 51, 151, 251, ...
    i = 100*(k-1) + 51
    nc_r[k, 1:6] = data_nc_r[i, 1:6]
    nc_r[k, 7] = r_pred[i]
}

p_nc_r = nc_r %>% 
    ggplot() +
    geom_line(aes(x = r, y = V, color = "red")) +
    geom_line(aes(x = r, y = V_pred)) +
    theme(legend.position = "none") + 
    ggtitle("r vs. V") +
    xlab("r, risk-free rate") +
    ylab("V, option price")


####################################
# non-constant: sigma, volatility
####################################

data_nc_sigma = data.frame("r" = rep(NA, n_obs_ep),
                       "sigma" = sigma_ep,
                       "D" = rep(NA, n_obs_ep),
                       "K" = rep(NA, n_obs_ep),
                       "S" = rep(NA, n_obs_ep),
                       "V" = rep(NA, n_obs_ep))


for (k in 1:n_sets_ep){
    # i:j goes over 100 indices per each k
    i = 100*(k-1) + 1
    j = 100*k
    data_nc_sigma$r[i:j] = rep(cp$r, 100)
    data_nc_sigma$sigma[i:j] = rep(sigma_ep[k], 100)
    data_nc_sigma$D[i:j] = rep(cp$D, 100)
    data_nc_sigma$K[i:j] = rep(cp$K, 100)
    data_nc_sigma$S[i:j] = seq(0, 3*cp$K, length.out = 100)
    data_nc_sigma$V[i:j] = fds_1s_A(time, cp$r, cp$K, sigma_ep[k], cp$D)
}

data_nc_sigma = as.matrix(data_nc_sigma)
sigma_pred = amer_mod %>% predict(data_nc_sigma[,1:5])
plot(sigma_pred, data_nc_sigma[,6], main = "sigma, non-constant", xlab = "actual (V)", ylab = "predicted (V)", cex = 0.3)
abline(0, 1, col = "red")

sigma_mae  = mean(abs((sigma_pred - data_nc_sigma[,6])))
sigma_mae
# 1.008805

# get median S for each set of 100
nc_sigma = data.frame("r" = rep(NA, n_sets_ep),
                  "sigma" = rep(NA, n_sets_ep),
                  "D" = rep(NA, n_sets_ep),
                  "K" = rep(NA, n_sets_ep),
                  "S" = rep(NA, n_sets_ep),
                  "V" = rep(NA, n_sets_ep),
                  "V_pred" = rep(NA, n_sets_ep))

for (k in 1:n_sets_ep){
    # k goes from 1 to 100
    # i = 51, 151, 251, ...
    i = 100*(k-1) + 51
    nc_sigma[k, 1:6] = data_nc_sigma[i, 1:6]
    nc_sigma[k, 7] = sigma_pred[i]
}

p_nc_sigma = nc_sigma %>% 
    ggplot() +
    geom_line(aes(x = sigma, y = V, color = "red")) +
    geom_line(aes(x = sigma, y = V_pred)) +
    theme(legend.position = "none") + 
    ggtitle("sigma vs. V") +
    xlab("sigma, volatility") +
    ylab("V, option price")

####################################
# non-constant: D, dividend yield
####################################
data_nc_D = data.frame("r" = rep(NA, n_obs_ep),
                       "sigma" = rep(NA, n_obs_ep),
                       "D" = D_ep,
                       "K" = rep(NA, n_obs_ep),
                       "S" = rep(NA, n_obs_ep),
                       "V" = rep(NA, n_obs_ep))


for (k in 1:n_sets_ep){
    # i:j goes over 100 indices per each k
    i = 100*(k-1) + 1
    j = 100*k
    data_nc_D$r[i:j] = rep(cp$r, 100)
    data_nc_D$sigma[i:j] = rep(cp$sigma, 100)
    data_nc_D$D[i:j] = rep(D_ep[k], 100)
    data_nc_D$K[i:j] = rep(cp$K, 100)
    data_nc_D$S[i:j] = seq(0, 3*cp$K, length.out = 100)
    data_nc_D$V[i:j] = fds_1s_A(time, cp$r, cp$K, cp$sigma, D_ep[k])
}

data_nc_D = as.matrix(data_nc_D)
D_pred = amer_mod %>% predict(data_nc_D[,1:5])
plot(D_pred, data_nc_D[,6], main = "D, non-constant", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

D_mae  = mean(abs((D_pred - data_nc_D[,6])))
D_mae
# 0.9009106 

# get median S for each set of 100
nc_D = data.frame("r" = rep(NA, n_sets_ep),
                      "sigma" = rep(NA, n_sets_ep),
                      "D" = rep(NA, n_sets_ep),
                      "K" = rep(NA, n_sets_ep),
                      "S" = rep(NA, n_sets_ep),
                      "V" = rep(NA, n_sets_ep),
                      "V_pred" = rep(NA, n_sets_ep))

for (k in 1:n_sets_ep){
    # k goes from 1 to 100
    # i = 51, 151, 251, ...
    i = 100*(k-1) + 51
    nc_D[k, 1:6] = data_nc_D[i, 1:6]
    nc_D[k, 7] = D_pred[i]
}

p_nc_D = nc_D %>% 
    ggplot() +
    geom_line(aes(x = D, y = V, color = "red")) +
    geom_line(aes(x = D, y = V_pred)) +
    theme(legend.position = "none") + 
    ggtitle("D vs. V") +
    xlab("D, dividend yield") +
    ylab("V, option price")

####################################
# non-constant: K, strike price
####################################
data_nc_K = data.frame("r" = rep(NA, n_obs_ep),
                       "sigma" = rep(NA, n_obs_ep),
                       "D" = rep(NA, n_obs_ep),
                       "K" = K_ep,
                       "S" = rep(NA, n_obs_ep),
                       "V" = rep(NA, n_obs_ep))


for (k in 1:n_sets_ep){
    # i:j goes over 100 indices per each k
    i = 100*(k-1) + 1
    j = 100*k
    data_nc_K$r[i:j] = rep(cp$r, 100)
    data_nc_K$sigma[i:j] = rep(cp$sigma, 100)
    data_nc_K$D[i:j] = rep(cp$D, 100)
    data_nc_K$K[i:j] = rep(K_ep[k], 100)
    data_nc_K$S[i:j] = seq(0, 3*K_ep[k], length.out = 100)
    data_nc_K$V[i:j] = fds_1s_A(time, cp$r, K_ep[k], cp$sigma, cp$D)
}

data_nc_K = as.matrix(data_nc_K)
K_pred = amer_mod %>% predict(data_nc_K[,1:5])
plot(K_pred, data_nc_K[,6], main = "K, non-constant", xlab = "predicted", ylab = "actual")
abline(0, 1, col = "red")

K_mae  = mean(abs((K_pred - data_nc_K[,6])))
K_mae
# 0.9870275

# get median S for each set of 100
nc_K = data.frame("r" = rep(NA, n_sets_ep),
                  "sigma" = rep(NA, n_sets_ep),
                  "D" = rep(NA, n_sets_ep),
                  "K" = rep(NA, n_sets_ep),
                  "S" = rep(NA, n_sets_ep),
                  "V" = rep(NA, n_sets_ep),
                  "V_pred" = rep(NA, n_sets_ep))

for (k in 1:n_sets_ep){
    # k goes from 1 to 100
    # i = 51, 151, 251, ...
    i = 100*(k-1) + 51
    nc_K[k, 1:6] = data_nc_K[i, 1:6]
    nc_K[k, 7] = K_pred[i]
}

p_nc_K = nc_K %>% 
    ggplot() +
    geom_line(aes(x = K, y = V, color = "red")) +
    geom_line(aes(x = K, y = V_pred)) +
    theme(legend.position = "none") + 
    ggtitle("K vs. V") +
    xlab("K, strike price") +
    ylab("V, option price")

################################################################################################
# II: Extrapolation Test Sets
################################################################################################

n_tst1 = 0.2 * obs_in_data
nn_tst1 = n_tst1 / 100

set.seed(42)
r_tst1 = runif(n_tst1, 0.035, 0.085)
sigma_tst1 = runif(n_tst1, 0.35, 0.85)
D_tst1 = runif(n_tst1, 0.035, 0.085)
K_tst1 = runif(n_tst1, 50, 100)

nstep_time1 = round(max(sigma_tst1), 1)^2 * 100^2   #6400
time1 = seq(0, 1, length.out = nstep_time1)

S_tst1 = c(rep(NA,n_tst1))
data_tst1 = data.frame("r" = rep(NA, n_tst1), 
                       "sigma" = rep(NA, n_tst1), 
                       "D" = rep(NA, n_tst1), 
                       "K" = rep(NA, n_tst1), 
                       "S" = rep(NA, n_tst1),
                       "V" = rep(NA, n_tst1))

for (k in 1:nn_tst1){
    # i:j goes over 100 indices per each k
    i = 100*(k-1) + 1
    j = 100*k
    print(k)
    print(i)
    print(j)
    data_tst1$r[i:j] = rep(r_tst1[k], 100)
    data_tst1$sigma[i:j] = rep(sigma_tst1[k], 100)
    data_tst1$D[i:j] = rep(D_tst1[k], 100)
    data_tst1$K[i:j] = rep(K_tst1[k], 100)
    data_tst1$S[i:j] = seq(0, 3*K_tst1[k], length.out = 100)
    data_tst1$V[i:j] = fds_1s_A(time1, r_tst1[k], K_tst1[k], sigma_tst1[k], D_tst1[k])
}

# mean(is.nan(data_tst1$V))

set.seed(42)
tst1_pred = amer_mod %>% predict(as.matrix(data_tst1[,1:5]))
tst1_act = as.matrix(data_tst1[,6])
plot(tst1_pred, tst1_act, main = "", xlab = "actual", ylab = "predicted")
abline(0, 1, col = "red")

tst1_rmse  = mean(abs((tst1_pred - tst1_act)))
tst1_rmse
# 2.574011