

setwd("C:/~")
rm(list=ls())

source("generateSample.R")
source("PenalizedLS.R")
source("functionalTS.R")
source("OracleTS.R")

###### local linear fit for functional observations ######
local.fit <- function(W, X.time){
  x.fit <- NULL
  for(j in 1 : length(W)){
    Wj <- W[[j]]
    wj.time <- X.time[[j]]
    tgrid.len <- length(wj.time)
    sample_size <- nrow(Wj)
    xj.fit <- matrix(NA, sample_size, tgrid.len)
    for (i in 1 : sample_size){
      wij <- as.vector(Wj[i, ])
      Wij.data <- data.frame(wij)
      Wij.data$t <- wj.time
      xij.fit <- locpol::locpol(wij ~ t, Wij.data, xeval = wj.time, bw = 0.05)
      xj.fit[i, ] <- xij.fit$lpFit$wij
    }
    x.fit[[j]] <- xj.fit
  }
  return(x.fit)
}

######## the codes for computing the PSR and NSR ###############
selection_Rate_f <- function(beta_true, beta_est, alpha_true, alpha_est){
  fun <- function(x) apply(x, 1, function(x) sum(abs(x)))
  indictor_true <- t(sapply(beta_true, fun))
  nonzero_location_true <- which(indictor_true != 0)
  zero_location_true <- which(indictor_true == 0)
  nonzero_true_num <- length(nonzero_location_true)
  zero_true_num <- length(zero_location_true)
  indictor_est <- t(sapply(beta_est, fun))
  nonzero_location_est <- which(indictor_est != 0)
  zero_location_est <- which(indictor_est == 0)
  Positive_num = length(base::intersect(nonzero_location_est, nonzero_location_true))
  PSR <- Positive_num/nonzero_true_num
  negative_num = length(base::intersect(zero_location_est, zero_location_true))
  NSR <- negative_num/zero_true_num
  return(list(PSR = PSR, NSR = NSR))
}

selection_Rate_s <- function(alpha_true, alpha_est){
  nonzero_location_true <- which(alpha_true != 0)
  zero_location_true <- which(alpha_true == 0)
  nonzero_true_num <- length(nonzero_location_true)
  zero_true_num <- length(zero_location_true)
  nonzero_location_est <- which(alpha_est != 0)
  zero_location_est <- which(alpha_est == 0)
  Positive_num = length(base::intersect(nonzero_location_est, nonzero_location_true))
  PSR <- Positive_num/nonzero_true_num
  negative_num = length(base::intersect(zero_location_est, zero_location_true))
  NSR <- negative_num/zero_true_num
  return(list(PSR = PSR, NSR = NSR))
}

######################################
########### simulation 1 #############

## generate AR sample and EX sample ##
AR_sample <- function(n_train, n_test, m, d, p, CorrNum, sigma2, rho){
  data_train <- ARSample(n_train, m, d, p, CorrNum, sigma2, rho)
  data_test <- ARSample(n_test, m, d, p, CorrNum, sigma2, rho)
  return(list(data.train = data_train, data.test = data_test))
}

EX_sample <- function(n_train, n_test, m, d, p, CorrNum, sigma2, rho){
  data_train <- EXSample(n_train, m, d, p, CorrNum, sigma2, rho)
  data_test <- EXSample(n_test, m, d, p, CorrNum, sigma2, rho)
  return(list(data.train = data_train, data.test = data_test))
}

## estimation with specified tuning parameters for given data ##
estimation_sn <- function(Y_train, W_train, Z_train, W_test, X_test,
                          Z_test, X.time, beta_true, alpha_true, sn){
  n_test <- nrow(Z_test)
  m <- ncol(Y_train)
  d <- length(W_train)
  trunPar <- rep(sn, d)
  Xhat_train <- local.fit(W_train, X.time)
  Xhat_test <- local.fit(W_test, X.time)
  pls <- PLSfun(Y_train, Xhat_train, Z_train, X.time, trunPar)
  fts <- FTSfun(Y_train, Xhat_train, Z_train, X.time, trunPar)
  beta_pls <- pls$beta
  beta_fts <- fts$beta
  alpha_pls <- pls$alpha
  alpha_fts <- fts$alpha
  SelectRate_f_pls <- selection_Rate_f(beta_true, beta_pls)
  SelectRate_f_fts <- selection_Rate_f(beta_true, beta_fts)
  SelectRate_s_pls <- selection_Rate_s(alpha_true, alpha_pls)
  SelectRate_s_fts <- selection_Rate_s(alpha_true, alpha_fts)
  PSR_f_pls <- SelectRate_f_pls$PSR
  NSR_f_pls <- SelectRate_f_pls$NSR
  PSR_f_fts <- SelectRate_f_fts$PSR
  NSR_f_fts <- SelectRate_f_fts$NSR
  PSR_s_pls <- SelectRate_s_pls$PSR
  NSR_s_pls <- SelectRate_s_pls$NSR
  PSR_s_fts <- SelectRate_s_fts$PSR
  NSR_s_fts <- SelectRate_s_fts$NSR
  MSEf_pls <- NULL
  MSEf_fts <- NULL
  fx_pls <- array(NA, c(n_test, m, d))
  fx_fts <- array(NA, c(n_test, m, d))
  fx_test <- array(NA, c(n_test, m, d))
  for (j in 1 : d){
    MSEf_pls[j] <- sum((beta_true[[j]] - beta_pls[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_pls[ , , j] <- Xhat_test[[j]] %*% t(beta_pls[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_fts[j] <- sum((beta_true[[j]] - beta_fts[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_fts[ , , j] <- Xhat_test[[j]] %*% t(beta_fts[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    fx_test[ , , j] <- X_test[[j]] %*% t(beta_true[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
  }
  Y_test <- apply(fx_test, c(1, 2), sum) + Z_test %*% alpha_true
  Y_pls <- apply(fx_pls, c(1, 2), sum) + Z_test %*% alpha_pls
  Y_fts <- apply(fx_fts, c(1, 2), sum) + Z_test %*% alpha_fts
  PE_pls <- mean((Y_pls - Y_test)^2) * m
  PE_fts <- mean((Y_fts - Y_test)^2) * m
  MSEf_pls <- sum(MSEf_pls)
  MSEf_fts <- sum(MSEf_fts)
  MSEs_pls <- sum((alpha_true - alpha_pls)^2)
  MSEs_fts <- sum((alpha_true - alpha_fts)^2)
  return(list(beta.pls = beta_pls, beta.fts = beta_fts, MSEf.pls = MSEf_pls,
              MSEf.fts = MSEf_fts, MSEs.pls = MSEs_pls, MSEs.fts = MSEs_fts, 
              PSRf.pls = PSR_f_pls, PSRf.fts = PSR_f_fts, NSRf.pls = NSR_f_pls, 
              NSRf.fts = NSR_f_fts, PE.pls = PE_pls, PE.fts = PE_fts, 
              alpha.pls = alpha_pls, alpha.fts = alpha_fts, PSRs.pls = PSR_s_pls,
              PSRs.fts = PSR_s_fts, NSRs.pls = NSR_s_pls, NSRs.fts = NSR_s_fts))
}

## estimation with tuning parameters selected by ABIC criterion for given data ##
estimation <- function(Y_train, W_train, Z_train, W_test, X_test,
                       Z_test, X.time, beta_true, alpha_true){
  n_test <- nrow(Z_test)
  m <- ncol(Y_train)
  d <- length(W_train)
  Xhat_train <- local.fit(W_train, X.time)
  Xhat_test <- local.fit(W_test, X.time)
  pls <- PLSfun.ABIC(Y_train, Xhat_train, Z_train, X.time)
  fts <- FTSfun.ABIC(Y_train, Xhat_train, Z_train, X.time)
  beta_pls <- pls$beta
  beta_fts <- fts$beta
  alpha_pls <- pls$alpha
  alpha_fts <- fts$alpha
  SelectRate_f_pls <- selection_Rate_f(beta_true, beta_pls)
  SelectRate_f_fts <- selection_Rate_f(beta_true, beta_fts)
  SelectRate_s_pls <- selection_Rate_s(alpha_true, alpha_pls)
  SelectRate_s_fts <- selection_Rate_s(alpha_true, alpha_fts)
  PSR_f_pls <- SelectRate_f_pls$PSR
  NSR_f_pls <- SelectRate_f_pls$NSR
  PSR_f_fts <- SelectRate_f_fts$PSR
  NSR_f_fts <- SelectRate_f_fts$NSR
  PSR_s_pls <- SelectRate_s_pls$PSR
  NSR_s_pls <- SelectRate_s_pls$NSR
  PSR_s_fts <- SelectRate_s_fts$PSR
  NSR_s_fts <- SelectRate_s_fts$NSR
  fx_test <- array(NA, c(n_test, m, d))
  fx_pls <- array(NA, c(n_test, m, d))
  fx_fts <- array(NA, c(n_test, m, d))
  MSEf_pls <- NULL
  MSEf_fts <- NULL
  for (j in 1 : d){
    MSEf_pls[j] <- sum((beta_true[[j]] - beta_pls[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_pls[ , , j] <- Xhat_test[[j]] %*% t(beta_pls[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_fts[j] <- sum((beta_true[[j]] - beta_fts[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_fts[ , , j] <- Xhat_test[[j]] %*% t(beta_fts[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    fx_test[ , , j] <- X_test[[j]] %*% t(beta_true[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
  }
  Y_test <- apply(fx_test, c(1, 2), sum) + Z_test %*% alpha_true
  Y_pls <- apply(fx_pls, c(1, 2), sum) + Z_test %*% alpha_pls
  Y_fts <- apply(fx_fts, c(1, 2), sum) + Z_test %*% alpha_fts
  PE_pls <- mean((Y_pls - Y_test)^2) * m
  PE_fts <- mean((Y_fts - Y_test)^2) * m
  MSEf_pls <- sum(MSEf_pls)
  MSEf_fts <- sum(MSEf_fts)
  MSEs_pls <- sum((alpha_true - alpha_pls)^2)
  MSEs_fts <- sum((alpha_true - alpha_fts)^2)
  return(list(beta.pls = beta_pls, beta.fts = beta_fts, MSEf.pls = MSEf_pls,
              MSEf.fts = MSEf_fts, MSEs.pls = MSEs_pls, snj.pls = pls$snj.aic,
              MSEs.fts = MSEs_fts, PSRf.pls = PSR_f_pls, PSRf.fts = PSR_f_fts,
              NSRf.pls = NSR_f_pls, NSRf.fts = NSR_f_fts, PE.pls = PE_pls,
              PE.fts = PE_fts,  snj.fts = fts$snj.aic,
              alpha.pls = alpha_pls, alpha.fts = alpha_fts, PSRs.pls = PSR_s_pls,
              PSRs.fts = PSR_s_fts, NSRs.pls = NSR_s_pls, NSRs.fts = NSR_s_fts))
}

#################################################################################
#set.seed(0)
Re = 200 # repeat 200 times
AR_data <- replicate(Re, AR_sample(n_train = 400, n_test = 500, m = 5, d = 4,
                                   p = 10, CorrNum = 3, sigma2 = 1, rho = 0.5))
#EX_data <- replicate(Re, EX_sample(n_train = 400, n_test = 500, m = 5, d = 4,
#                                   p = 10, CorrNum = 3, sigma2 = 1, rho = 0.5))
#save(AR_data, file = "./AR_data.RData")

sn_range <- c(1 : 8, 10, 12, 15, 16, 18, 20)
R = length(sn_range)
result_pls <- matrix(NA, 14, 11)
result_fts <- matrix(NA, 14, 11)
for (j in 1 : length(sn_range)){
  MSEf_pls <- rep(NA, Re)
  MSEf_fts <- rep(NA, Re)
  MSEs_pls <- rep(NA, Re)
  MSEs_fts <- rep(NA, Re)
  PE_pls <- rep(NA, Re)
  PE_fts <- rep(NA, Re)
  PSRf_pls <- rep(NA, Re)
  PSRf_fts <- rep(NA, Re)
  NSRf_pls <- rep(NA, Re)
  NSRf_fts <- rep(NA, Re)
  PSRs_pls <- rep(NA, Re)
  PSRs_fts <- rep(NA, Re)
  NSRs_pls <- rep(NA, Re)
  NSRs_fts <- rep(NA, Re)
  snj_pls_list <- list()
  snj_fts_list <- list()
  sn <- sn_range[j]
  for (i in 1 : Re){
    data_train <- AR_data[, i]$data.train
    data_test <- AR_data[, i]$data.test
    Y_train <- data_train$Y
    X_train <- data_train$X
    Z_train <- data_train$Z
    W_train <- data_train$W
    X.time <- data_train$X.time
    beta_true <- data_train$beta
    alpha_true <- data_train$alpha
    X_test <- data_test$X
    Z_test <- data_test$Z
    W_test <- data_test$W
    #AR <- estimation_sn(Y_train, W_train, Z_train, W_test, X_test,
    #                    Z_test, X.time, beta_true, alpha_true, sn)
    AR <- estimation(Y_train, W_train, Z_train, W_test, X_test,
                     Z_test, X.time, beta_true, alpha_true)
    MSEf_pls[i] <- AR$MSEf.pls
    MSEf_fts[i] <- AR$MSEf.fts
    MSEs_pls[i] <- AR$MSEs.pls
    MSEs_fts[i] <- AR$MSEs.fts
    PE_pls[i] <- AR$PE.pls
    PE_fts[i] <- AR$PE.fts
    PSRf_pls[i] <- AR$PSRf.pls
    PSRf_fts[i] <- AR$PSRf.fts
    NSRf_pls[i] <- AR$NSRf.pls
    NSRf_fts[i] <- AR$NSRf.fts
    PSRs_pls[i] <- AR$PSRs.pls
    PSRs_fts[i] <- AR$PSRs.fts
    NSRs_pls[i] <- AR$NSRs.pls
    NSRs_fts[i] <- AR$NSRs.fts
    snj_pls_list[[i]] <- AR$snj.pls
    snj_fts_list[[i]] <- AR$snj.fts
  }
  result_pls[j, ] <- c(sn, mean(PSRf_pls), mean(NSRf_pls), mean(PSRs_pls), mean(NSRs_pls),
                       mean(MSEf_pls), sd(MSEf_pls), mean(MSEs_pls), sd(MSEs_pls),
                       mean(PE_pls), sd(PE_pls))
  result_fts[j, ] <- c(sn, mean(PSRf_fts), mean(NSRf_fts), mean(PSRs_fts), mean(NSRs_fts),
                       mean(MSEf_fts), sd(MSEf_fts), mean(MSEs_fts), sd(MSEs_fts),
                       mean(PE_fts), sd(PE_fts))
}

colnames(result_pls) <- c("sn", "PSRf", "NSRf", "PSRs", "NSRs", "MSEf_mean", "MSEf_sd",
                          "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")
colnames(result_fts) <- c("sn", "PSRf", "NSRf", "PSRs", "NSRs", "MSEf_mean", "MSEf_sd",
                          "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")

snj.fts <- unlist(snj_fts_list)
snj.fts.mean <- sum(snj.fts)/sum(snj.fts != 0)
snj.pls <- unlist(snj_pls_list)
snj.pls.mean <- sum(snj.pls)/sum(snj.pls != 0)
result_pls <- rbind(result_pls, c(snj.pls.mean, mean(PSRf_pls), mean(NSRf_pls), mean(PSRs_pls), mean(NSRs_pls),
                                  mean(MSEf_pls), sd(MSEf_pls), mean(MSEs_pls), sd(MSEs_pls),
                                  mean(PE_pls), sd(PE_pls)))
result_fts <- rbind(result_fts, c(snj.fts.mean, mean(PSRf_fts), mean(NSRf_fts), mean(PSRs_fts), mean(NSRs_fts),
                                  mean(MSEf_fts), sd(MSEf_fts), mean(MSEs_fts), sd(MSEs_fts),
                                  mean(PE_fts), sd(PE_fts)))

write.csv(result_fts, "ARresult_fts_n400sn.csv")
write.csv(result_pls, "ARresult_pls_n400sn.csv")


######################################
######### simulation 2 #########

## simulation for AR samples with tuning parameters selected by ABIC procedure ##
AR_simulation <- function(n_train, n_test, m, d, p, CorrNum, sigma2, rho){
  data_train <- ARSample(n_train, m, d, p, CorrNum, sigma2, rho)
  beta_true <- data_train$beta
  alpha_true <- data_train$alpha
  Y_train <- data_train$Y
  Z_train <- data_train$Z
  W_train <- data_train$W
  X.time <- data_train$X.time
  Xhat_train <- local.fit(W_train, X.time)
  pls <- PLSfun.ABIC(Y_train, Xhat_train, Z_train, X.time)
  fts <- FTSfun.ABIC(Y_train, Xhat_train, Z_train, X.time)
  Orts <- OrTS.AIC(Y_train, Xhat_train, Z_train, X.time, beta_true, alpha_true)
  Orls <- OrLS.AIC(Y_train, Xhat_train, Z_train, X.time, beta_true, alpha_true)
  beta_pls <- pls$beta
  beta_fts <- fts$beta
  beta_Orts <- Orts$beta
  beta_Orls <- Orls$beta
  alpha_pls <- pls$alpha
  alpha_fts <- fts$alpha
  alpha_Orts <- Orts$alpha
  alpha_Orls <- Orls$alpha
  SelectRate_f_pls <- selection_Rate_f(beta_true, beta_pls)
  SelectRate_f_fts <- selection_Rate_f(beta_true, beta_fts)
  SelectRate_s_pls <- selection_Rate_s(alpha_true, alpha_pls)
  SelectRate_s_fts <- selection_Rate_s(alpha_true, alpha_fts)
  PSR_f_pls <- SelectRate_f_pls$PSR
  NSR_f_pls <- SelectRate_f_pls$NSR
  PSR_f_fts <- SelectRate_f_fts$PSR
  NSR_f_fts <- SelectRate_f_fts$NSR
  PSR_s_pls <- SelectRate_s_pls$PSR
  NSR_s_pls <- SelectRate_s_pls$NSR
  PSR_s_fts <- SelectRate_s_fts$PSR
  NSR_s_fts <- SelectRate_s_fts$NSR
  data_test <- ARSample(n_test, m, d, p, CorrNum, sigma2, rho)
  X_test <- data_test$X
  Z_test <- data_test$Z
  W_test <- data_test$W
  Xhat_test <- local.fit(W_test, X.time)
  fx_test <- array(NA, c(n_test, m, d))
  fx_pls <- array(NA, c(n_test, m, d))
  fx_fts <- array(NA, c(n_test, m, d))
  fx_Orts <- array(NA, c(n_test, m, d))
  fx_Orls <- array(NA, c(n_test, m, d))
  MSEf_pls <- NULL
  MSEf_fts <- NULL
  MSEf_Orts <- NULL
  MSEf_Orls <- NULL
  for (j in 1 : d){
    MSEf_pls[j] <- sum((beta_true[[j]] - beta_pls[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_pls[ , , j] <- Xhat_test[[j]] %*% t(beta_pls[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_fts[j] <- sum((beta_true[[j]] - beta_fts[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_fts[ , , j] <- Xhat_test[[j]] %*% t(beta_fts[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_Orts[j] <- sum((beta_true[[j]] - beta_Orts[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_Orts[ , , j] <- Xhat_test[[j]] %*% t(beta_Orts[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_Orls[j] <- sum((beta_true[[j]] - beta_Orls[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_Orls[ , , j] <- Xhat_test[[j]] %*% t(beta_Orls[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    fx_test[ , , j] <- X_test[[j]] %*% t(beta_true[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
  }
  Y_test <- apply(fx_test, c(1, 2), sum) + Z_test %*% alpha_true
  Y_pls <- apply(fx_pls, c(1, 2), sum) + Z_test %*% alpha_pls
  Y_fts <- apply(fx_fts, c(1, 2), sum) + Z_test %*% alpha_fts
  Y_Orts <- apply(fx_Orts, c(1, 2), sum) + Z_test %*% alpha_Orts
  Y_Orls <- apply(fx_Orls, c(1, 2), sum) + Z_test %*% alpha_Orls
  PE_pls <- mean((Y_pls - Y_test)^2) * m
  PE_fts <- mean((Y_fts - Y_test)^2) * m
  PE_Orts <- mean((Y_Orts - Y_test)^2) * m
  PE_Orls <- mean((Y_Orls - Y_test)^2) * m
  MSEf_pls <- sum(MSEf_pls)
  MSEf_fts <- sum(MSEf_fts)
  MSEf_Orts <- sum(MSEf_Orts)
  MSEf_Orls <- sum(MSEf_Orls)
  MSEs_pls <- sum((alpha_true - alpha_pls)^2)
  MSEs_fts <- sum((alpha_true - alpha_fts)^2)
  MSEs_Orts <- sum((alpha_true - alpha_Orts)^2)
  MSEs_Orls <- sum((alpha_true - alpha_Orls)^2)
  return(list(beta.pls = beta_pls,  beta.fts = beta_fts, beta.Orts = beta_Orts,
              alpha.pls = alpha_pls, alpha.fts = alpha_fts, alpha.Orts = alpha_Orts,
              MSEf.pls = MSEf_pls, MSEf.fts = MSEf_fts, MSEf.Orts = MSEf_Orts, 
              MSEs.pls = MSEs_pls, MSEs.fts = MSEs_fts, MSEs.Orts = MSEs_Orts, 
              PSRf.pls = PSR_f_pls, PSRf.fts = PSR_f_fts, 
              NSRf.pls = NSR_f_pls, NSRf.fts = NSR_f_fts,
              PSRs.pls = PSR_s_pls, PSRs.fts = PSR_s_fts, 
              NSRs.pls = NSR_s_pls, NSRs.fts = NSR_s_fts,
              PE.pls = PE_pls, PE.fts = PE_fts, PE.Orts = PE_Orts,
              MSEf.Orls = MSEf_Orls, MSEs.Orls = MSEs_Orls, PE.Orls = PE_Orls))
}

result_pls <- matrix(NA, 10, 11)
result_fts <- matrix(NA, 10, 11)
result_Orts <- matrix(NA, 10, 7)
result_Orls <- matrix(NA, 10, 7)
MSEf_list <- list()
MSEs_list <- list()
PE_list <- list()
R <- 200
rho_vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
for (j in 1 : 10){
  AR1 <- replicate(R, AR_simulation(n_train = 400, n_test = 500, m = 5, d = 4,
                                    p = 10, CorrNum = 3, sigma2 = 1, rho = rho_vec[j]))
  MSEf_pls <- rep(NA, R)
  MSEf_fts <- rep(NA, R)
  MSEf_Orts <- rep(NA, R)
  MSEf_Orls <- rep(NA, R)
  MSEs_pls <- rep(NA, R)
  MSEs_fts <- rep(NA, R)
  MSEs_Orts <- rep(NA, R)
  MSEs_Orls <- rep(NA, R)
  PE_pls <- rep(NA, R)
  PE_fts <- rep(NA, R)
  PE_Orts <- rep(NA, R)
  PE_Orls <- rep(NA, R)
  PSRf_pls <- rep(NA, R)
  PSRf_fts <- rep(NA, R)
  NSRf_pls <- rep(NA, R)
  NSRf_fts <- rep(NA, R)
  PSRs_pls <- rep(NA, R)
  PSRs_fts <- rep(NA, R)
  NSRs_pls <- rep(NA, R)
  NSRs_fts <- rep(NA, R)
  for (i in 1 : R){
    MSEf_pls[i] <- AR1[, i]$MSEf.pls
    MSEf_fts[i] <- AR1[, i]$MSEf.fts
    MSEf_Orts[i] <- AR1[, i]$MSEf.Orts
    MSEf_Orls[i] <- AR1[, i]$MSEf.Orls
    MSEs_pls[i] <- AR1[, i]$MSEs.pls
    MSEs_fts[i] <- AR1[, i]$MSEs.fts
    MSEs_Orts[i] <- AR1[, i]$MSEs.Orts
    MSEs_Orls[i] <- AR1[, i]$MSEs.Orls
    PE_pls[i] <- AR1[, i]$PE.pls
    PE_fts[i] <- AR1[, i]$PE.fts
    PE_Orts[i] <- AR1[, i]$PE.Orts
    PE_Orls[i] <- AR1[, i]$PE.Orls
    PSRf_pls[i] <- AR1[, i]$PSRf.pls
    PSRf_fts[i] <- AR1[, i]$PSRf.fts
    NSRf_pls[i] <- AR1[, i]$NSRf.pls
    NSRf_fts[i] <- AR1[, i]$NSRf.fts
    PSRs_pls[i] <- AR1[, i]$PSRs.pls
    PSRs_fts[i] <- AR1[, i]$PSRs.fts
    NSRs_pls[i] <- AR1[, i]$NSRs.pls
    NSRs_fts[i] <- AR1[, i]$NSRs.fts
  }
  result_pls[j, ] <- c(rho_vec[j], mean(PSRf_pls), mean(NSRf_pls), mean(PSRs_pls), mean(NSRs_pls),
                       mean(MSEf_pls), sd(MSEf_pls), mean(MSEs_pls), sd(MSEs_pls),
                       mean(PE_pls), sd(PE_pls))
  result_fts[j, ] <- c(rho_vec[j], mean(PSRf_fts), mean(NSRf_fts), mean(PSRs_fts), mean(NSRs_fts),
                       mean(MSEf_fts), sd(MSEf_fts), mean(MSEs_fts), sd(MSEs_fts),
                       mean(PE_fts), sd(PE_fts))
  result_Orts[j, ] <- c(rho_vec[j], mean(MSEf_Orts), sd(MSEf_Orts), mean(MSEs_Orts), sd(MSEs_Orts),
                        mean(PE_Orts), sd(PE_Orts))
  result_Orls[j, ] <- c(rho_vec[j], mean(MSEf_Orls), sd(MSEf_Orls), mean(MSEs_Orls), sd(MSEs_Orls),
                        mean(PE_Orls), sd(PE_Orls))
  MSEf_list[[j]] <- cbind(MSEf_pls, MSEf_fts, MSEf_Orls, MSEf_Orts)
  MSEs_list[[j]] <- cbind(MSEs_pls,  MSEs_fts, MSEs_Orls, MSEs_Orts)
  PE_list[[j]] <- cbind(PE_pls, PE_fts, PE_Orls, PE_Orts)
  print(j)
}

colnames(result_pls) <- c("rho", "PSRf", "NSRf", "PSRs", "NSRs", "MSEf_mean", "MSEf_sd",
                          "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")
colnames(result_fts) <- c("rho", "PSRf", "NSRf", "PSRs", "NSRs", "MSEf_mean", "MSEf_sd",
                          "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")
colnames(result_Orts) <- c("rho", "MSEf_mean", "MSEf_sd", "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")
colnames(result_Orls) <- c("rho", "MSEf_mean", "MSEf_sd", "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")

write.csv(result_fts, "AR_result_fts_n400.csv", quote = FALSE)
write.csv(result_pls, "AR_result_pls_n400.csv", quote = FALSE)
write.csv(result_Orls, "AR_result_Orls_n400.csv", quote = FALSE)
write.csv(result_Orts, "AR_result_Orts_n400.csv", quote = FALSE)

colname1 <- rep(rho_vec, each = 800)
colname2 <- rep(c("PLS", "fTS", "LS_oracle", "TS_oracle"), each =200)
colname2 <- rep(colname2, length(rho_vec))
MSEf_data <- unlist(MSEf_list)
MSEf_data <- cbind(colname1, colname2, MSEf_data)
colnames(MSEf_data) <- c("rho", "methods", "MSEf_value")
write.csv(MSEf_data, "AR_MSEfn400.csv", quote = FALSE, row.names = FALSE)

MSEs_data <- unlist(MSEs_list)
MSEs_data <- cbind(colname1, colname2, MSEs_data)
colnames(MSEs_data) <- c("rho", "methods", "MSEs_value")
write.csv(MSEs_data, "AR_MSEsn400.csv", quote = FALSE, row.names = FALSE)

PE_data <- unlist(PE_list)
PE_data <- cbind(colname1, colname2, PE_data)
colnames(PE_data) <- c("rho", "methods", "PE_value")
write.csv(PE_data, "AR_PEn400.csv", quote = FALSE, row.names = FALSE)


## simulation for EX samples with tuning parameters selected by ABIC procedure ##
EX_simulation <- function(n_train, n_test, m, d, p, CorrNum, sigma2, rho){
  data_train <- EXSample(n_train, m, d, p, CorrNum, sigma2, rho)
  beta_true <- data_train$beta
  alpha_true <- data_train$alpha
  Y_train <- data_train$Y
  Z_train <- data_train$Z
  W_train <- data_train$W
  X.time <- data_train$X.time
  Xhat_train <- local.fit(W_train, X.time)
  pls <- PLSfun.ABIC(Y_train, Xhat_train, Z_train, X.time)
  fts <- FTSfun.ABIC(Y_train, Xhat_train, Z_train, X.time)
  Orts <- OrTS.AIC(Y_train, Xhat_train, Z_train, X.time, beta_true, alpha_true)
  Orls <- OrLS.AIC(Y_train, Xhat_train, Z_train, X.time, beta_true, alpha_true)
  beta_pls <- pls$beta
  beta_fts <- fts$beta
  beta_Orts <- Orts$beta
  beta_Orls <- Orls$beta
  alpha_pls <- pls$alpha
  alpha_fts <- fts$alpha
  alpha_Orts <- Orts$alpha
  alpha_Orls <- Orls$alpha
  SelectRate_f_pls <- selection_Rate_f(beta_true, beta_pls)
  SelectRate_f_fts <- selection_Rate_f(beta_true, beta_fts)
  SelectRate_s_pls <- selection_Rate_s(alpha_true, alpha_pls)
  SelectRate_s_fts <- selection_Rate_s(alpha_true, alpha_fts)
  PSR_f_pls <- SelectRate_f_pls$PSR
  NSR_f_pls <- SelectRate_f_pls$NSR
  PSR_f_fts <- SelectRate_f_fts$PSR
  NSR_f_fts <- SelectRate_f_fts$NSR
  PSR_s_pls <- SelectRate_s_pls$PSR
  NSR_s_pls <- SelectRate_s_pls$NSR
  PSR_s_fts <- SelectRate_s_fts$PSR
  NSR_s_fts <- SelectRate_s_fts$NSR
  data_test <- EXSample(n_test, m, d, p, CorrNum, sigma2, rho)
  X_test <- data_test$X
  Z_test <- data_test$Z
  W_test <- data_test$W
  Xhat_test <- local.fit(W_test, X.time)
  fx_test <- array(NA, c(n_test, m, d))
  fx_pls <- array(NA, c(n_test, m, d))
  fx_fts <- array(NA, c(n_test, m, d))
  fx_Orts <- array(NA, c(n_test, m, d))
  fx_Orls <- array(NA, c(n_test, m, d))
  MSEf_pls <- NULL
  MSEf_fts <- NULL
  MSEf_Orts <- NULL
  MSEf_Orls <- NULL
  for (j in 1 : d){
    MSEf_pls[j] <- sum((beta_true[[j]] - beta_pls[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_pls[ , , j] <- Xhat_test[[j]] %*% t(beta_pls[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_fts[j] <- sum((beta_true[[j]] - beta_fts[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_fts[ , , j] <- Xhat_test[[j]] %*% t(beta_fts[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_Orts[j] <- sum((beta_true[[j]] - beta_Orts[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_Orts[ , , j] <- Xhat_test[[j]] %*% t(beta_Orts[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    MSEf_Orls[j] <- sum((beta_true[[j]] - beta_Orls[[j]])^2) * (X.time[[j]][2] - X.time[[j]][1])
    fx_Orls[ , , j] <- Xhat_test[[j]] %*% t(beta_Orls[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
    fx_test[ , , j] <- X_test[[j]] %*% t(beta_true[[j]]) * (X.time[[j]][2] - X.time[[j]][1])
  }
  Y_test <- apply(fx_test, c(1, 2), sum) + Z_test %*% alpha_true
  Y_pls <- apply(fx_pls, c(1, 2), sum) + Z_test %*% alpha_pls
  Y_fts <- apply(fx_fts, c(1, 2), sum) + Z_test %*% alpha_fts
  Y_Orts <- apply(fx_Orts, c(1, 2), sum) + Z_test %*% alpha_Orts
  Y_Orls <- apply(fx_Orls, c(1, 2), sum) + Z_test %*% alpha_Orls
  PE_pls <- mean((Y_pls - Y_test)^2) * m
  PE_fts <- mean((Y_fts - Y_test)^2) * m
  PE_Orts <- mean((Y_Orts - Y_test)^2) * m
  PE_Orls <- mean((Y_Orls - Y_test)^2) * m
  MSEf_pls <- sum(MSEf_pls)
  MSEf_fts <- sum(MSEf_fts)
  MSEf_Orts <- sum(MSEf_Orts)
  MSEf_Orls <- sum(MSEf_Orls)
  MSEs_pls <- sum((alpha_true - alpha_pls)^2)
  MSEs_fts <- sum((alpha_true - alpha_fts)^2)
  MSEs_Orts <- sum((alpha_true - alpha_Orts)^2)
  MSEs_Orls <- sum((alpha_true - alpha_Orls)^2)
  return(list(beta.pls = beta_pls, beta.fts = beta_fts, beta.Orts = beta_Orts,
              alpha.pls = alpha_pls, alpha.fts = alpha_fts, alpha.Orts = alpha_Orts,
              MSEf.pls = MSEf_pls, MSEf.fts = MSEf_fts, MSEf.Orts = MSEf_Orts, 
              MSEf.Orls = MSEf_Orls, MSEs.pls = MSEs_pls, MSEs.fts = MSEs_fts, 
              MSEs.Orts = MSEs_Orts, MSEs.Orls = MSEs_Orls, PSRf.pls = PSR_f_pls, 
              PSRf.fts = PSR_f_fts, NSRf.pls = NSR_f_pls, NSRf.fts = NSR_f_fts,
              PSRs.pls = PSR_s_pls, PSRs.fts = PSR_s_fts, NSRs.pls = NSR_s_pls, NSRs.fts = NSR_s_fts,
              PE.pls = PE_pls, PE.fts = PE_fts, PE.Orts = PE_Orts, PE.Orls = PE_Orls))
}

#set.seed(22)
result_pls <- matrix(NA, 10, 11)
result_fts <- matrix(NA, 10, 11)
result_Orts <- matrix(NA, 10, 7)
result_Orls <- matrix(NA, 10, 7)
MSEf_list <- list()
MSEs_list <- list()
PE_list <- list()
R <- 200
rho_vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
for (j in 1 : 10){
  EX1 <- replicate(R, EX_simulation(n_train = 200, n_test = 500, m = 5, d = 4,
                                    p = 10, CorrNum = 3, sigma2 = 1, rho = rho_vec[j]))
  MSEf_pls <- rep(NA, R)
  MSEf_fts <- rep(NA, R)
  MSEf_Orts <- rep(NA, R)
  MSEf_Orls <- rep(NA, R)
  MSEs_pls <- rep(NA, R)
  MSEs_fts <- rep(NA, R)
  MSEs_Orts <- rep(NA, R)
  MSEs_Orls <- rep(NA, R)
  PE_pls <- rep(NA, R)
  PE_fts <- rep(NA, R)
  PE_Orts <- rep(NA, R)
  PE_Orls <- rep(NA, R)
  PSRf_pls <- rep(NA, R)
  PSRf_fts <- rep(NA, R)
  NSRf_pls <- rep(NA, R)
  NSRf_fts <- rep(NA, R)
  PSRs_pls <- rep(NA, R)
  PSRs_fts <- rep(NA, R)
  NSRs_pls <- rep(NA, R)
  NSRs_fts <- rep(NA, R)
  for (i in 1 : R){
    MSEf_pls[i] <- EX1[, i]$MSEf.pls
    MSEf_fts[i] <- EX1[, i]$MSEf.fts
    MSEf_Orts[i] <- EX1[, i]$MSEf.Orts
    MSEf_Orls[i] <- EX1[, i]$MSEf.Orls
    MSEs_pls[i] <- EX1[, i]$MSEs.pls
    MSEs_fts[i] <- EX1[, i]$MSEs.fts
    MSEs_Orts[i] <- EX1[, i]$MSEs.Orts
    MSEs_Orls[i] <- EX1[, i]$MSEs.Orls
    PE_pls[i] <- EX1[, i]$PE.pls
    PE_fts[i] <- EX1[, i]$PE.fts
    PE_Orts[i] <- EX1[, i]$PE.Orts
    PE_Orls[i] <- EX1[, i]$PE.Orls
    PSRf_pls[i] <- EX1[, i]$PSRf.pls
    PSRf_fts[i] <- EX1[, i]$PSRf.fts
    NSRf_pls[i] <- EX1[, i]$NSRf.pls
    NSRf_fts[i] <- EX1[, i]$NSRf.fts
    PSRs_pls[i] <- EX1[, i]$PSRs.pls
    PSRs_fts[i] <- EX1[, i]$PSRs.fts
    NSRs_pls[i] <- EX1[, i]$NSRs.pls
    NSRs_fts[i] <- EX1[, i]$NSRs.fts
  }
  result_pls[j, ] <- c(rho_vec[j], mean(PSRf_pls), mean(NSRf_pls), mean(PSRs_pls), mean(NSRs_pls), 
                       mean(MSEf_pls), sd(MSEf_pls), mean(MSEs_pls), sd(MSEs_pls), 
                       mean(PE_pls), sd(PE_pls))
  result_fts[j, ] <- c(rho_vec[j], mean(PSRf_fts), mean(NSRf_fts), mean(PSRs_fts), mean(NSRs_fts), 
                       mean(MSEf_fts), sd(MSEf_fts), mean(MSEs_fts), sd(MSEs_fts), 
                       mean(PE_fts), sd(PE_fts))
  result_Orts[j, ] <- c(rho_vec[j], mean(MSEf_Orts), sd(MSEf_Orts), mean(MSEs_Orts), sd(MSEs_Orts), 
                       mean(PE_Orts), sd(PE_Orts))
  result_Orls[j, ] <- c(rho_vec[j], mean(MSEf_Orls), sd(MSEf_Orls), mean(MSEs_Orls), sd(MSEs_Orls), 
                       mean(PE_Orls), sd(PE_Orls))
  MSEf_list[[j]] <- cbind(MSEf_pls, MSEf_fts, MSEf_Orls, MSEf_Orts)
  MSEs_list[[j]] <- cbind(MSEs_pls, MSEs_fts, MSEs_Orls, MSEs_Orts)
  PE_list[[j]] <- cbind(PE_pls, PE_fts, PE_Orls, PE_Orts)
  print(j)
}

colnames(result_pls) <- c("rho", "PSRf", "NSRf", "PSRs", "NSRs", "MSEf_mean", "MSEf_sd", 
                          "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")
colnames(result_fts) <- c("rho", "PSRf", "NSRf", "PSRs", "NSRs", "MSEf_mean", "MSEf_sd", 
                          "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")
colnames(result_Orts) <- c("rho", "MSEf_mean", "MSEf_sd", "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")
colnames(result_Orls) <- c("rho", "MSEf_mean", "MSEf_sd", "MSEs_mean", "MSEs_sd", "PE_mean", "PE_sd")

write.csv(result_fts, "EX_result_fts_n200.csv", quote = FALSE)
write.csv(result_pls, "EX_result_pls_n200.csv", quote = FALSE)
write.csv(result_Orls, "EX_result_Orls_n200.csv", quote = FALSE)
write.csv(result_Orts, "EX_result_Orts_n200.csv", quote = FALSE)

colname1 <- rep(rho_vec, each = 800)
colname2 <- rep(c("PLS", "fTS", "LS_oracle", "TS_oracle"), each =200)
colname2 <- rep(colname2, length(rho_vec))
MSEf_data <- unlist(MSEf_list)
MSEf_data <- cbind(colname1, colname2, MSEf_data)
colnames(MSEf_data) <- c("rho", "methods", "MSEf_value")
write.csv(MSEf_data, "EX_MSEfn200.csv", quote = FALSE, row.names = FALSE)

MSEs_data <- unlist(MSEs_list)
MSEs_data <- cbind(colname1, colname2, MSEs_data)
colnames(MSEs_data) <- c("rho", "methods", "MSEs_value")
write.csv(MSEs_data, "EX_MSEsn200.csv", quote = FALSE, row.names = FALSE)

PE_data <- unlist(PE_list)
PE_data <- cbind(colname1, colname2, PE_data)
colnames(PE_data) <- c("rho", "methods", "PE_value")
write.csv(PE_data, "EX_PEn200.csv", quote = FALSE, row.names = FALSE)

