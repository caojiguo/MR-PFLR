
## initial values of the truncation parameters for functional predictors ##
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

###### Two-stage algorithm for estimating beta and Omega #######
TSalgo <- function(Y, XZ, n, m, X_group, Z_group){ 
  y <- as.vector(t(Y))   
  group <- c(X_group, Z_group)
  ## stage I ##
  betafit1 <- grpreg::grpreg(XZ, y, penalty = "grSCAD", group = group, lambda.min = 1e-12)
  coef1 <- as.vector(grpreg::select(betafit1, "BIC")$beta)[-1]
  Yfit1 <- t(matrix(XZ %*% coef1, m, n))
  S1 <- cor(Yfit1 - Y)
  Omegafit1 <- GGMncv::ggmncv(S1, n, penalty = "scad", progress = FALSE)
  Omega1 <- Omegafit1$Theta
  ## stage II ##
  Omega1_root <- pracma::rootm(Omega1, p = 2)$B
  Lambda_root <- kronecker(diag(1, n), Omega1_root)
  y_tilde <- Lambda_root %*% y
  XZ_tilde <- Lambda_root %*% XZ
  betafit2 <- grpreg::grpreg(XZ_tilde, y_tilde, penalty = "grSCAD", group = group, lambda.min = 1e-12)
  coef2 <- as.vector(grpreg::select(betafit2, "BIC")$beta)[-1]
  # betafit2 <- grpreg::cv.grpreg(XZ_tilde, y_tilde, penalty = "grSCAD", group = group_new)
  # coef2 <- as.vector(coef(betafit2))[-1]
  Yfit2 <- t(matrix(XZ %*% coef2, m, n))
  # Yfit2 <- t(matrix(betafit2$y, m, n))
  # S2 <- cor(Yfit2 - Y)
  # Omegafit2 <- GGMncv::ggmncv(S2, n, penalty = "scad", progress = FALSE, lambda_min_ratio = 1e-4)
  # Omega2 <- Omegafit2$Theta
  beta <- coef2[1 : length(X_group)]
  alpha <- coef2[-c(1 : length(X_group))]
  return(list(betaPar = beta, alpha = alpha, fit = Yfit2, fit0 = Yfit1))
} 


################################################################
## functional two-stage procedure for estimating beta with specified truncation parameters ##
FTSfun <- function(Y, X, Z, X.tgrid, trunPar){
  if(is.list(X)) X = X else X = list(X)
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  d <- length(X)
  p <- ncol(Z)
  y <- as.vector(t(Y))
  snjVec <- trunPar  # the truncation parameter vector for d functional predictors
  snj_sum <- sum(snjVec)
  repNum <- rep(snjVec, m)
  X_group <- rep(1 : (m * d), repNum)
  Z_group <- (m * d + 1) : ((m * d) + m * p)
  group <- c(X_group, Z_group)
  x.basis <- NULL
  Scores_est <- NULL
  for (j in 1 : d){
    xfdat <- fda.usc::fdata(X[[j]], X.tgrid[[j]])
    xpc <- fda.usc::fdata2pc(xfdat, ncomp = snjVec[j])
    x.basis[[j]] <- xpc$rotation$data # the estimated fpc basis
    Scores_est[[j]] <- xpc$x[, 1 : snjVec[j]] # the first snj estimated fpc scores
  }
  Scores_est <- matrix(unlist(Scores_est), ncol = snj_sum) 
  krofun <- function(x) kronecker(diag(1, m), x)
  X1 <- t(matrix(apply(Scores_est, 1, krofun), m * dim(Scores_est)[2],
                m * dim(Scores_est)[1])) # design matrix of fpc scores 
  Z1 <- t(matrix(apply(Z, 1, krofun), nrow = p * m)) # design matrix of scalar covariates
  XZ <- cbind(X1, Z1)
  TSfit <- TSalgo(Y, XZ, n, m, X_group, Z_group)
  betaPar_est<- TSfit$betaPar
  betaPar_est_list <- split(betaPar_est, X_group)
  beta_est <- NULL
  for (j in 1 : d){
    beta_est_j <- matrix(NA, m, length(X.tgrid[[j]]))
    for (l in 1 : m){
      beta_est_j[l, ] <- betaPar_est_list[[(l-1)*d+j]] %*% x.basis[[j]] 
    }
    beta_est[[j]] <- beta_est_j
  }
  alpha_est <- matrix(TSfit$alpha, nrow = p)
  Yfit <- TSfit$fit
  Yfit0 <- TSfit$fit0
  Omega <- TSfit$Omega
  return(list(beta = beta_est, alpha = alpha_est, fit = Yfit, fit0 = Yfit0, Omega = Omega))
}

### fTS-based AIC criterion for selecting s_nj ###
FTSfun.AIC <- function(Y, X, Z_select, n, m, d, X.tgrid, snj.max, beta.bic){
  #snj.max: the maximal truncation parameters for all j
  #beta.bic: a list of all estimated beta based on FTS_bic procedure 
  y <- as.vector(t(Y))
  snjVec <- rep(snj.max, d)
  snj_sum <- sum(snjVec) 
  beta_bic <- beta.bic
  fun <- function(x) apply(x, 1, function(x) sum(abs(x)))
  indictor <- t(sapply(beta_bic, fun)) # d * m matrix with 0 indicating location of zero beta
  nonzero_location <- which(indictor != 0)
  nonzero_indictor <- which(indictor != 0, arr.ind = TRUE) #col corresponds to m locations
  nonzero_j_all <- unique(nonzero_indictor[, 1]) # all indices of nonzero beta in the direction of j  
  nonzero_l_all <- unique(nonzero_indictor[, 2]) # all indices of nonzero beta in the direction of l 
  nonzero_p_num <- length(nonzero_j_all)
  nonzero_m_num <- length(nonzero_l_all)
  if(any(nonzero_p_num == 0 || nonzero_m_num == 0) ) stop('No coefficients exist')
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
  X_group <- rep(1 : (m * d), repNum)
  W_list <- lapply(split(W, X_group), matrix, nrow = n * m) # a list consists of d*m matrices
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
    W_new_list[[s]] <- W_select
    snjMa_list[[s]] <- snj_Ma
    XZ_select <- cbind(W_select, Z_select)
    ## stage I ##
    glsfit1 <- lm(y ~ XZ_select - 1)
    betaPar1 <- glsfit1$coefficients
    Yfit1 <- t(matrix(glsfit1$fitted.values, m, n))
    S1 <- cor(Yfit1 - Y)
    Omegafit1 <- GGMncv::ggmncv(S1, n, penalty = "scad", progress = FALSE)
    Omega1 <- Omegafit1$Theta
    ## stage II ##
    Omega1_root <- pracma::rootm(Omega1, p = 2)$B
    Lambda_root <- kronecker(diag(1, n), Omega1_root)
    y_tilde <- Lambda_root %*% y
    XZ_tilde <- Lambda_root %*% XZ_select
    glsfit2 <- lm(y_tilde ~ XZ_tilde - 1)
    betaPar2 <- glsfit2$coefficients
    Yfit2 <- t(matrix(XZ_select %*% betaPar2, m, n))
    rss <- sum((Yfit2 - Y)^2)
    AIC_scores[s] <- log(rss) + 2 * sum(snj_Ma) /n
  }
  s_location <- which.min(AIC_scores)
  snj_est <- snjMa_list[[s_location]] # truncation parameter in (j,l)th location
  W_new <- W_new_list[[s_location]]
  return(list(snj = snj_est, W_new = W_new))
}

## FPWLS method for estimation and variable selection with 
##  truncation parameters selected by ABIC criterion ##
FTSfun.ABIC <- function(Y, X, Z, X.tgrid){
  if(is.list(X)) X = X else X = list(X)
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  d <- length(X)
  p <- ncol(Z)
  snj_ini <- iniTrun(X, X.tgrid)  # initial truncation parameter vector for d functional predictors
  FTSfit.bic <- FTSfun(Y, X, Z, X.tgrid, snj_ini)
  beta_bic <- FTSfit.bic$beta
  alpha_bic <- as.vector(FTSfit.bic$alpha)
  alp_nonzero_location <- which(alpha_bic != 0)
  krofun <- function(x) kronecker(diag(1, m), x)
  Z1 <- t(matrix(apply(Z, 1, krofun), nrow = p * m))
  Z_select <- Z1[, alp_nonzero_location]
  snj_max <- max(snj_ini)
  FTSfit_aic <- FTSfun.AIC(Y, X, Z_select, n, m, d, X.tgrid, snj_max, beta_bic)
  snj <- FTSfit_aic$snj
  nonzero_location <- which(snj != 0)
  snj_select <- as.vector(snj)[nonzero_location]
  X_group_new <- rep(1 : sum(snj != 0), snj_select)
  Z_group_new <- (sum(snj != 0) + 1) : (sum(snj != 0) + length(alp_nonzero_location))
  W_new <- FTSfit_aic$W_new
  XZ_new <- cbind(W_new, Z_select)
  ftsfit <- TSalgo(Y, XZ_new, n, m, X_group_new, Z_group_new)
  alpha_est <- rep(0, p * m)
  alpha_est[alp_nonzero_location] <- ftsfit$alpha
  alpha_est <- matrix(alpha_est, nrow = p)
  betaPar_est <- ftsfit$betaPar
  betaPar_est_list <- split(betaPar_est, X_group_new)
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
  coef <- c(betaPar_est, ftsfit$alpha)
  Yfit <- t(matrix(XZ_new %*% coef, m, n))
  return(list(beta = beta_est, alpha = alpha_est, fit = Yfit, snj.ini = snj_ini, snj.aic = snj))
}