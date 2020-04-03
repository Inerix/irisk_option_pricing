#####################
# @ Params
# sig: volatilty
# cor: Correlation between underlying
# r: Risk free interest rate
# K: Strike
# S: stock price
# w: weight
# weights: Vector representing weights of the stocks in the basket
#####################

# data setup
data_2s = read.csv("data/two_stock.csv")
colnames(data_2s) = c("sig1", "sig2", "cor", "r", "K", "w1", "w2", "S1", "S2", "V")
n_2s = nrow(data_2s)
trn_2s_idx = sample(seq(from = 1, to = n_2s, by = 1), 0.8*n_2s)
data_2s_trn = as.matrix(data_2s[trn_2s_idx,-])
data_2s_tst = as.matrix(data_2s[-trn_2s_idx,])

basket_mod = keras_model_sequential() %>% 
    layer_dense(units = 64, activation = "relu", input_shape = 5) %>% 
    layer_dense(units = 1) %>% 
    
    compile(
        optimizer = "rmsprop",
        loss = "mse",
        metrics = c("mae")  
    )

t3 = proc.time()
set.seed(42)
amer_fit = amer_mod %>%
    fit(
        x = data_trn[, 1:5],
        y = data_trn[, 6],
        epochs = 500, 
        batch = 1024,  
        validation_split = 0.25,
        verbose = TRUE
    )
t4 = proc.time()
t4-t3



