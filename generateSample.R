
X.phi <- function(NumberOfphi, tgrid){
  phi <- matrix(NA, NumberOfphi, length(tgrid))
  for (j in 1 : NumberOfphi)
    if (j %% 2 == 0){
      phi[j, ] <- 2 * j^(-1) * sqrt(2) * sin((j - 1) * pi * tgrid)
    } else {
      phi[j, ] <- 2 * j^(-1) * sqrt(2) * cos(j * pi * tgrid)
    }
  return(phi)
}

# AR covariance matrix
ar_cor <- function(m, sigma2, rho) {
  exponent <- abs(matrix(1:m - 1, nrow = m, ncol = m, byrow = TRUE) - 
                    (1:m - 1))
  sigma2 * rho^exponent
}
# EX covariance matrix
ex_cor <- function(m, sigma2, rho) {
  matrix(sigma2 * rho, nrow = m, ncol = m) + diag(sigma2 - sigma2 * rho, m) 
}

xz.cov <- function(r, p){
  xz.cov <- matrix(NA, 4, p)
  for (i in 1 : 4){
    for (j in 1 : p){
      xz.cov[i, j] <- r^(abs(i - j) + 1)
    }
  }
  return(xz.cov)
}

Vfun <- function(NumberOfphi, tgrid){
  xsi <- rnorm(NumberOfphi, 0, 1)
  V <- colSums(xsi * X.phi(NumberOfphi, tgrid)) 
  return(V)
}

Zfun <- function(n, NumberOfphi, tgrid, r, p){
  data <- list()
  cov.xi <- diag(1, 4)
  cov.z <- ar_cor(p, 1, 0.5)
  cov.xz <- xz.cov(r, p)
  cov1 <- cbind(cov.xi, cov.xz)
  cov2 <- cbind(t(cov.xz), cov.z)
  covXZ <- rbind(cov1, cov2)
  xsi4.Zq <- mvtnorm::rmvnorm(n, mean = rep(0, p + 4), sigma = covXZ)
  xsi1 <- xsi4.Zq[, 1 : 4]
  xsi2 <- mvtnorm::rmvnorm(n, mean = rep(0, NumberOfphi - 4), sigma = diag(1, NumberOfphi-4)) 
  xsi <- cbind(xsi1, xsi2)
  V2 <- xsi %*% X.phi(NumberOfphi, tgrid)
  data$Z <- xsi4.Zq[, -(1 : 4)]
  data$V2 <- V2
  return(data)
}

beta.Par <- function(m, NumberOfphi, d){
  betaPar <- array(0, c(m, NumberOfphi, d))
  for (l in 1 : 2){
    betaPar[l, 1 : 4 , 1] <- c(1, -0.8, 0.6, 1.2)
  }
  for (l in 3 : 4){
    for (k in 1 : 3){
    betaPar[l, k, 2] <- (-1)^(k + 1) * (1.2 - 0.2 * k)
    }
    for (k in 4 : NumberOfphi){
      betaPar[l, k, 2] <- 8 * (-1)^(k) * (k - 1)^(-4)
    }
  }
  return(betaPar)
}

alpha.fun <- function(p, m){
  alpha = matrix(0, p, m)
  alpha[, 1] <- c(1, 1, rep(0, p - 2))
  alpha[, 2] <- c(rep(0, p - 1), 1.4)
  alpha[, 4] <- c(1.2, 0, 1.5, rep(0, p - 3))
  alpha[, 5] <- c(0, 2, 0.8, rep(0, p - 3))
  return(alpha)
}

#################################################################
############# generate data for AR error ###########
ARSample <- function(n, m, d, p, CorrNum, sigma2, rho){
  tgrid <- seq(0, 1, len = 100)
  NumberOfphi = 50
  data <- list()
  v.list <- NULL
  for (j in 1 : (d + 2)){
    v.list[[j]] <- t(replicate(n, Vfun(NumberOfphi, tgrid)))
  }
  ZV2 <- Zfun(n, NumberOfphi, tgrid, r = 0.2, p)
  v.list[[3]] <- ZV2$V2
  x <- NULL
  w <- NULL
  x.time <- NULL
  for (j in 2 : (CorrNum + 1)){  #CorrNum: number of dependent predictors
    x[[j-1]] <- v.list[[j]] + 0.5 * (v.list[[j-1]] + v.list[[j+1]])
    x.error <- rnorm(n * length(tgrid), 0, 1)
    w[[j-1]] <- x[[j-1]] + matrix(x.error, nrow = n)
    x.time[[j-1]] <- tgrid
  }
  for (j in (CorrNum + 1) : d) {
    x[[j]] <- v.list[[j + 2]]
    x.error <- rnorm(n * length(tgrid), 0, 1)
    w[[j]] <- x[[j]] + matrix(x.error, nrow = n)
    x.time[[j]] <- tgrid
  }
  data$X <- x
  data$W <- w
  data$X.time <- x.time
  data$Z <- ZV2$Z
  betaPar <- beta.Par(m, NumberOfphi, d)
  beta <- NULL
  fx.array <- array(NA, c(n, m, d))
  for (j in 1 : d){
    beta[[j]] <- betaPar[ , , j] %*% X.phi(NumberOfphi, x.time[[j]])    
    fx.array[ , , j] <- x[[j]] %*% t(beta[[j]]) * (x.time[[j]][2] - x.time[[j]][1])
  }
  data$beta <- beta
  alpha <- alpha.fun(p, m)
  data$alpha <- alpha
  errorCov <- ar_cor(m, sigma2, rho)
  eps <- mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = errorCov)
  y <- apply(fx.array, c(1, 2), sum) + (ZV2$Z) %*% alpha + eps
  data$Y <- y
  data$errorCov <- errorCov
  return(data)
}

#######################################################
############# generate data for EX error ###########
EXSample <- function(n, m, d, p, CorrNum, sigma2, rho){
  tgrid <- seq(0, 1, len = 100)
  NumberOfphi = 50
  data <- list()
  v.list <- NULL
  for (j in 1 : (d + 2)){
    v.list[[j]] <- t(replicate(n, Vfun(NumberOfphi, tgrid)))
  }
  ZV2 <- Zfun(n, NumberOfphi, tgrid, r=0.2, p)
  v.list[[3]] <- ZV2$V2
  x <- NULL
  w <- NULL
  x.time <- NULL
  for (j in 2 : (CorrNum + 1)){  #CorrNum: number of dependent predictors
    x[[j-1]] <- v.list[[j]] + 0.5 * (v.list[[j-1]] + v.list[[j+1]])
    x.error <- rnorm(n * length(tgrid), 0, 1)
    w[[j-1]] <- x[[j-1]] + matrix(x.error, nrow = n)
    x.time[[j-1]] <- tgrid
  }
  for (j in (CorrNum + 1) : d) {
    x[[j]] <- v.list[[j + 2]]
    x.error <- rnorm(n * length(tgrid), 0, 1)
    w[[j]] <- x[[j]] + matrix(x.error, nrow = n)
    x.time[[j]] <- tgrid
  }
  data$X <- x
  data$W <- w
  data$X.time <- x.time
  data$Z <- ZV2$Z
  betaPar <- beta.Par(m, NumberOfphi, d)
  beta <- NULL
  fx.array <- array(NA, c(n, m, d))
  for (j in 1 : d){
    beta[[j]] <- betaPar[ , , j] %*% X.phi(NumberOfphi, x.time[[j]])    
    fx.array[ , , j] <- x[[j]] %*% t(beta[[j]]) * (x.time[[j]][2] - x.time[[j]][1])
  }
  data$beta <- beta
  alpha <- alpha.fun(p, m)
  data$alpha <- alpha
  errorCov <- ex_cor(m, sigma2, rho)
  eps <- mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = errorCov)
  y <- apply(fx.array, c(1, 2), sum) + (ZV2$Z) %*% alpha + eps
  data$Y <- y
  data$errorCov <- errorCov
  return(data)
}