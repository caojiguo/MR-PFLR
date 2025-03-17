
##### initial values of the truncation parameters for predictors #####
iniTrun <- function(x.fit, x.tgrid){
  d <- length(x.fit)
  trunPar <- rep(1, d)
  for (j in 1 : d){
    xfdat <- fda.usc::fdata(x.fit[[j]], x.tgrid[[j]])
    xpc <- fda.usc::fdata2pc(xfdat)
    x.sd <- xpc$d
    x.var <- x.sd^2/sum(x.sd^2)
    x.cumvar <- cumsum(x.var)
    trunPar[j] <- which(x.cumvar >= 0.99)[1]
  }
  return(trunPar)
}

########################################################################
## FPLS-based oracle estimation method with truncation parameters selected by AIC criterion 
OrLS.AIC <- function(Y, X, Z, X.tgrid, beta0, alpha0){
  if(is.list(X)) X = X else X = list(X)
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  d <- length(X)
  p <- ncol(Z)
  y <- as.vector(t(Y))
  snj_ini <- iniTrun(X, X.tgrid)
  snj_max <- max(snj_ini)
  snjVec <- rep(snj_max, d)
  snj_sum <- sum(snjVec) 
  fun <- function(x) apply(x, 1, function(x) sum(abs(x)))
  beta_indictor <- t(sapply(beta0, fun)) # d * m matrix with 0 indicating location of zero beta
  nonzero_beta_location <- which(beta_indictor != 0)
  nonzero_beta_indictor <- which(beta_indictor != 0, arr.ind = TRUE) #col corresponds to m locations
  nonzero_beta_num <- length(nonzero_beta_location)
  nonzero_j_all <- unique(nonzero_beta_indictor[, 1]) # all indices of nonzero beta in the direction of j  
  nonzero_l_all <- unique(nonzero_beta_indictor[, 2]) # all indices of nonzero beta in the direction of l 
  nonzero_p_num <- length(nonzero_j_all)
  nonzero_m_num <- length(nonzero_l_all)
  nonzero_j_comb <- combn(unlist(sapply(snjVec[nonzero_j_all], function(x) 2 : x)),
                          nonzero_p_num) # all combinations of selected s_nj in nonzero predictors
  nonzero_j_comb <- nonzero_j_comb[ , !duplicated(t(nonzero_j_comb))] # remove duplicate groups
  xbasis <- list()
  Scores_estlist <- list()
  for (j in 1 : d){
    xfdat <- fda.usc::fdata(X[[j]], X.tgrid[[j]])
    xpc <- fda.usc::fdata2pc(xfdat, ncomp = snjVec[j])
    xbasis[[j]] <- xpc$rotation$data # the estimated fpc basis
    Scores_estlist[[j]] <- xpc$x[, 1 : snjVec[j]] # the first snj estimated fpc scores
  }
  Scores_est <- matrix(unlist(Scores_estlist), ncol = snj_sum) 
  krofun <- function(x) kronecker(diag(1, m), x)
  W <- t(matrix(apply(Scores_est, 1, krofun), m * dim(Scores_est)[2],
                m * dim(Scores_est)[1])) # design matrix of fpc scores 
  repNum <- rep((snjVec) * n * m, m)
  W_group <- rep(1 : (m * d), repNum)
  W_list <- lapply(split(W, W_group), matrix, nrow = n * m) # a list consists of p*m matrices
  Z1 <- t(matrix(apply(Z, 1, krofun), nrow = p * m)) # design matrix of scalar covariates
  nonzero_alpha_location <- which(alpha0 != 0)
  Z1_significant <- Z1[, nonzero_alpha_location]
  AIC_scores <- rep(NA, ncol(nonzero_j_comb))
  snjMa_list <- list()
  W_new_list <- list()
  coef_list <- list()
  for (s in 1 : ncol(nonzero_j_comb)){
    W_select_list <- list()
    snj_Ma <- matrix(0, nrow = d, ncol = m) # d * m 
    for (l in 1 : m){
      for (j in 1 : d){
        location_full <- (l-1)*d+j
        if (location_full %in% nonzero_beta_location) {
          location_j <- which(nonzero_j_all == j)
          snj_Ma[j, l] <- nonzero_j_comb[location_j, s]
          W_select_list[[location_full]] <- W_list[[location_full]][, 1 : nonzero_j_comb[location_j, s]]
        }
      }
    }
    W_select <- matrix(unlist(W_select_list), nrow = n * m)
    snjMa_list[[s]] <- snj_Ma
    WZ_select <- cbind(W_select, Z1_significant)
    lsfit <- lm(y ~ WZ_select - 1)
    coefPar <- lsfit$coefficients
    coef_list[[s]] <- coefPar
    Yfit <- t(matrix(WZ_select %*% coefPar, m, n))
    rss <- sum((Yfit - Y)^2)
    AIC_scores[s] <- log(rss) + 2 * sum(snj_Ma) / n
  }
  s_location <- which.min(AIC_scores)
  snj_est <- snjMa_list[[s_location]] # d*m matrix of truncation parameters
  coefPar <- coef_list[[s_location]]
  betaPar_len <- 1 : sum(snj_est)
  betaPar_est <- coefPar[betaPar_len]
  betaPar_group <- rep(1 : nonzero_beta_num, snj_est[snj_est != 0])
  betaPar_est_list <- split(betaPar_est, betaPar_group)
  beta_est <- NULL
  for (j in 1 : d){
    beta_est_j <- matrix(NA, m, length(X.tgrid[[j]]))
    for (l in 1 : m){
      location_full <- (l-1)*d+j
      if (location_full %in% nonzero_beta_location) {
        location_betaPar <- which(nonzero_beta_location == location_full)
        xfdat <- fda.usc::fdata(X[[j]], X.tgrid[[j]])
        xpc <- fda.usc::fdata2pc(xfdat, ncomp = snj_est[j, l])
        xbasis <- xpc$rotation$data 
        beta_est_j[l, ] <- betaPar_est_list[[location_betaPar]] %*% xbasis 
      }
      else beta_est_j[l, ] <- 0
    }
    beta_est[[j]] <- beta_est_j
  }
  alpha_est <- rep(0, p * m)
  alpha_est[nonzero_alpha_location] <- coefPar[-betaPar_len]
  alpha_est <- matrix(alpha_est, nrow = p)
  return(list(beta = beta_est, alpha = alpha_est))
}

## FPWLS-based oracle estimation method with truncation parameters selected by AIC criterion 
OrTS.AIC <- function(Y, X, Z, X.tgrid, beta0, alpha0){
  if(is.list(X)) X = X else X = list(X)
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  d <- length(X)
  p <- ncol(Z)
  y <- as.vector(t(Y))
  snj_ini <- iniTrun(X, X.tgrid)
  snj_max <- max(snj_ini)
  snjVec <- rep(snj_max, d)
  snj_sum <- sum(snjVec) 
  fun <- function(x) apply(x, 1, function(x) sum(abs(x)))
  beta_indictor <- t(sapply(beta0, fun)) # d * m matrix with 0 indicating location of zero beta
  nonzero_beta_location <- which(beta_indictor != 0)
  nonzero_beta_indictor <- which(beta_indictor != 0, arr.ind = TRUE) #col corresponds to m locations
  nonzero_beta_num <- length(nonzero_beta_location)
  nonzero_j_all <- unique(nonzero_beta_indictor[, 1]) # all indices of nonzero beta in the direction of j  
  nonzero_l_all <- unique(nonzero_beta_indictor[, 2]) # all indices of nonzero beta in the direction of l 
  nonzero_p_num <- length(nonzero_j_all)
  nonzero_m_num <- length(nonzero_l_all)
  nonzero_j_comb <- combn(unlist(sapply(snjVec[nonzero_j_all], function(x) 2 : x)),
                          nonzero_p_num) # all combinations of selected s_nj in nonzero predictors
  nonzero_j_comb <- nonzero_j_comb[ , !duplicated(t(nonzero_j_comb))] # remove duplicate groups
  xbasis <- list()
  Scores_estlist <- list()
  for (j in 1 : d){
    xfdat <- fda.usc::fdata(X[[j]], X.tgrid[[j]])
    xpc <- fda.usc::fdata2pc(xfdat, ncomp = snjVec[j])
    xbasis[[j]] <- xpc$rotation$data # the estimated fpc basis
    Scores_estlist[[j]] <- xpc$x[, 1 : snjVec[j]] # the first snj estimated fpc scores
  }
  Scores_est <- matrix(unlist(Scores_estlist), ncol = snj_sum) 
  krofun <- function(x) kronecker(diag(1, m), x)
  W <- t(matrix(apply(Scores_est, 1, krofun), m * dim(Scores_est)[2],
                m * dim(Scores_est)[1])) # design matrix of fpc scores 
  repNum <- rep((snjVec) * n * m, m)
  W_group <- rep(1 : (m * d), repNum)
  W_list <- lapply(split(W, W_group), matrix, nrow = n * m) # a list consists of d*m matrices
  Z1 <- t(matrix(apply(Z, 1, krofun), nrow = p * m)) # design matrix of scalar covariates
  nonzero_alpha_location <- which(alpha0 != 0)
  Z1_significant <- Z1[, nonzero_alpha_location]
  AIC_scores <- rep(NA, ncol(nonzero_j_comb))
  snjMa_list <- list()
  W_new_list <- list()
  coef_list <- list()
  for (s in 1 : ncol(nonzero_j_comb)){
    W_select_list <- list()
    snj_Ma <- matrix(0, nrow = d, ncol = m) # d * m 
    for (l in 1 : m){
      for (j in 1 : d){
        location_full <- (l-1)*d+j
        if (location_full %in% nonzero_beta_location) {
          location_j <- which(nonzero_j_all == j)
          snj_Ma[j, l] <- nonzero_j_comb[location_j, s]
          W_select_list[[location_full]] <- W_list[[location_full]][, 1 : nonzero_j_comb[location_j, s]]
        }
      }
    }
    W_select <- matrix(unlist(W_select_list), nrow = n * m)
    snjMa_list[[s]] <- snj_Ma
    WZ_select <- cbind(W_select, Z1_significant)
    ## stage I ##
    glsfit1 <- lm(y ~ WZ_select - 1)
    Yfit1 <- t(matrix(glsfit1$fitted.values, m, n))
    S1 <- cor(Yfit1 - Y)
    Omegafit1 <- GGMncv::ggmncv(S1, n, penalty = "scad", progress = FALSE)
    Omega1 <- Omegafit1$Theta
    ## stage II ##
    Omega1_root <- pracma::sqrtm(Omega1)$B
    Lambda_root <- kronecker(diag(1, n), Omega1_root)
    y_tilde <- Lambda_root %*% y
    WZ_tilde <- Lambda_root %*% WZ_select
    glsfit2 <- lm(y_tilde ~ WZ_tilde - 1)
    coefPar2 <- glsfit2$coefficients
    coef_list[[s]] <- coefPar2
    Yfit2 <- t(matrix(WZ_select %*% coefPar2, m, n))
    rss <- sum((Yfit2 - Y)^2)
    AIC_scores[s] <- log(rss) + 2 * sum(snj_Ma) / n
  }
  s_location <- which.min(AIC_scores)
  snj_est <- snjMa_list[[s_location]] # d*m matrix of truncation parameters
  coefPar <- coef_list[[s_location]]
  betaPar_len <- 1 : sum(snj_est)
  betaPar_est <- coefPar[betaPar_len]
  betaPar_group <- rep(1 : nonzero_beta_num, snj_est[snj_est != 0])
  betaPar_est_list <- split(betaPar_est, betaPar_group)
  beta_est <- NULL
  for (j in 1 : d){
    beta_est_j <- matrix(NA, m, length(X.tgrid[[j]]))
    for (l in 1 : m){ 
      location_full <- (l-1)*d+j
      if (location_full %in% nonzero_beta_location) {
        location_betaPar <- which(nonzero_beta_location == location_full)
        xfdat <- fda.usc::fdata(X[[j]], X.tgrid[[j]])
        xpc <- fda.usc::fdata2pc(xfdat, ncomp = snj_est[j, l])
        xbasis <- xpc$rotation$data 
        beta_est_j[l, ] <- betaPar_est_list[[location_betaPar]] %*% xbasis 
      }
      else beta_est_j[l, ] <- 0
    }
    beta_est[[j]] <- beta_est_j
  }
  alpha_est <- rep(0, p * m)
  alpha_est[nonzero_alpha_location] <- coefPar[-betaPar_len]
  alpha_est <- matrix(alpha_est, nrow = p)
  return(list(beta = beta_est, alpha = alpha_est))
}
