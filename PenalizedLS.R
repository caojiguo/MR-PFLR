
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
## PLS method for estimation and variable selection with specified truncation parameters 
PLSfun <- function(Y, X, Z, X.tgrid, trunPar){
  if(is.list(X)) X = X else X = list(X)
  if(is.matrix(Y)) {(n = nrow(Y)) && (m = ncol(Y))} else {(n = length(Y)) && (m = 1)}
  d <- length(X)
  p <- ncol(Z)
  y <- as.vector(t(Y))
  snjVec <- trunPar  # the truncation parameter vector for d functional predictors
  snj_sum <- sum(snjVec)
  repNum <- rep(snjVec, m)
  W_group <- rep(1 : (m * d), repNum)
  Z_group <- (m * d + 1) : ((m * d) + m * p)
  group <- c(W_group, Z_group)
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
  Z1 <- t(matrix(apply(Z, 1, krofun), nrow = p * m)) # design matrix of scalar covariates
  WZ <- cbind(W, Z1)
  BICfit <- grpreg::grpreg(WZ, y, penalty = "grSCAD", group = group, 
                            lambda.min = 1e-12)
  coef_est <- as.vector(grpreg::select(BICfit, "BIC")$beta)[-1]
  betaPar_est <- coef_est[1 : length(W_group)] 
  alpha_vec <- coef_est[-c(1 : length(W_group))]
  alpha_ma <- matrix(alpha_vec, nrow = p)  #with dimension p * m
  betaPar_estlist <- split(betaPar_est, W_group) 
  beta_est <- NULL
  for (j in 1 : d){
    beta_est_j <- matrix(NA, m, length(X.tgrid[[j]]))
    for (l in 1 : m){ 
      beta_est_j[l, ] <- betaPar_estlist[[(l-1)*d+j]] %*% xbasis[[j]] 
    }
    beta_est[[j]] <- beta_est_j
  }
  Yfit = t(matrix(WZ %*% coef_est, m, n))
  return(list(beta = beta_est, alpha = alpha_ma, fit = Yfit))
}

### FPLS-based AIC criterion for selecting s_nj ###
PLSfun.AIC <- function(Y, X, Z_select, n, m, d, X.tgrid, snj.max, beta.bic){
  #snj.max: the maximal truncation parameters for all j
  #beta.bic: a list of all estimated beta based on PLS_bic procedure 
  y <- as.vector(t(Y))
  snjVec <- rep(snj.max, d)
  snj_sum <- sum(snjVec) 
  beta_bic <- beta.bic
  fun <- function(x) apply(x, 1, function(x) sum(abs(x)))
  if(m == 1) {beta.indictor <- as.matrix(sapply(beta_bic, fun))} else {beta.indictor <- t(sapply(beta_bic, fun))}  # d * m matrix with 0 indicating location of zero beta
  nonzero_location <- which(beta.indictor != 0)
  nonzero_indictor <- which(beta.indictor != 0, arr.ind = TRUE) #col corresponds to m locations
  nonzero_j_all <- unique(nonzero_indictor[, 1]) # all indices of nonzero beta in the direction of j  
  nonzero_l_all <- unique(nonzero_indictor[, 2]) # all indices of nonzero beta in the direction of l 
  nonzero_p_num <- length(nonzero_j_all)
  nonzero_m_num <- length(nonzero_l_all)
  if(any(nonzero_p_num == 0 || nonzero_m_num == 0) ) stop('No coefficients exist')
  nonzero_j_comb <- combn(unlist(sapply(snjVec[nonzero_j_all], function(x) 2 : x)),
                          nonzero_p_num) # all combinations of selected s_nj in nonzero predictors
  nonzero_j_comb <- nonzero_j_comb[ , !duplicated(t(nonzero_j_comb))] # remove duplicate groups
  if(is.matrix(nonzero_j_comb)) {nonzero_j_comb = nonzero_j_comb} else {nonzero_j_comb = t(as.matrix(nonzero_j_comb))}
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
  W_list <- lapply(split(W, W_group), matrix, nrow = n * m)
  AIC_scores <- rep(NA, ncol(nonzero_j_comb))
  snjMa_list <- list()
  W_new_list <- list()
  for (s in 1 : ncol(nonzero_j_comb)){
    W_select_list <- list()
    snj_Ma <- matrix(0, nrow = d, ncol = m) # d * m 
    for (l in 1 : m){
      for (j in 1 : d){
        location_full <- (l-1)*d+j
        if (location_full %in% nonzero_location) {
          location_j <- which(nonzero_j_all == j)
          snj_Ma[j, l] <- nonzero_j_comb[location_j, s]
          W_select_list[[location_full]] <- W_list[[location_full]][, 1 : nonzero_j_comb[location_j, s]]
        }
      }
    }
    W_select <- matrix(unlist(W_select_list), nrow = n * m)
    snjMa_list[[s]] <- snj_Ma
    W_new_list[[s]] <- W_select
    WZ_select <- cbind(W_select, Z_select)
    olsfit <- lm(y ~ WZ_select - 1)
    res <- olsfit$residuals
    rss <- sum(res^2)
    AIC_scores[s] <- log(rss) + 2 * sum(snj_Ma) / n
  }
  s_location <- which.min(AIC_scores)
  snj_est <- snjMa_list[[s_location]] # truncation parameter in (j,l)th location
  W_new <- W_new_list[[s_location]]
  return(list(W.new = W_new, snj = snj_est))
}

## FPLS method for estimation and variable selection with 
##  truncation parameters selected by ABIC criterion ##
PLSfun.ABIC <- function(Y, X, Z, X.tgrid){
  if(is.list(X)) X = X else X = list(X)
  if(is.matrix(Y)) {(n = nrow(Y)) && (m = ncol(Y))} else {(n = length(Y)) && (m = 1)}
  d <- length(X)
  p <- ncol(Z)
  y <- as.vector(t(Y))
  snj_ini <- iniTrun(X, X.tgrid)  # initial truncation parameter vector for d functional predictors
  PLSfit_bic <-  PLSfun(Y, X, Z, X.tgrid, snj_ini)
  beta_bic <- PLSfit_bic$beta
  alpha_bic <- as.vector(PLSfit_bic$alpha)
  alp_nonzero_location <- which(alpha_bic != 0)
  krofun <- function(x) kronecker(diag(1, m), x)
  Z1 <- t(matrix(apply(Z, 1, krofun), nrow = p * m))
  Z_select <- Z1[, alp_nonzero_location]
  snj_max <- max(snj_ini)
  PLSfit_aic <-  PLSfun.AIC(Y, X, Z_select, n, m, d, X.tgrid, snj_max, beta_bic)
  snj <- PLSfit_aic$snj
  nonzero_location <- which(snj != 0)
  snj_select <- as.vector(snj)[nonzero_location]
  W_group <- rep(1 : sum(snj != 0), snj_select)
  if (sum(alp_nonzero_location) == 0){
    group_new <- W_group
    WZ_new <- PLSfit_aic$W.new
  }else{
    Z_group <- (sum(snj != 0) + 1) : (sum(snj != 0) + length(alp_nonzero_location))
    group_new <- c(W_group, Z_group)
    W_new <- PLSfit_aic$W.new
    WZ_new <- cbind(W_new, Z_select)
  }
  betafit <- grpreg::grpreg(WZ_new, y, penalty = "grSCAD", group = group_new, 
                            lambda.min = 1e-12)
  coef_est <- as.vector(grpreg::select(betafit, "BIC")$beta)[-1]
  alpha_est <- rep(0, p * m)
  alpha_est[alp_nonzero_location] <- coef_est[-c(1 : length(W_group))]
  alpha_est <- matrix(alpha_est, nrow = p)
  betaPar_est_list <- split(coef_est[1 : length(W_group)], W_group)
  beta_est <- NULL
  for (j in 1 : d){
    beta_est_j <- matrix(NA, m, length(X.tgrid[[j]]))
    for (l in 1 : m){ 
      location_full <- (l-1)*d+j
      if (location_full %in% nonzero_location) {
        location_betaPar <- which(nonzero_location == location_full)
        xfdat <- fda.usc::fdata(X[[j]], X.tgrid[[j]])
        xpc <- fda.usc::fdata2pc(xfdat, ncomp = snj[j, l])
        xbasis <- xpc$rotation$data 
        beta_est_j[l, ] <- betaPar_est_list[[location_betaPar]] %*% xbasis 
      }
      else beta_est_j[l, ] <- 0
    }
    beta_est[[j]] <- beta_est_j
  }
  Yfit <- t(matrix(WZ_new %*% coef_est, m, n))
  return(list(beta = beta_est, alpha = alpha_est, 
              fit = Yfit, snj.ini = snj_ini, snj.aic = snj))
}
