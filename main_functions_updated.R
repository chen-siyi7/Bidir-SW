library(MASS)  

first_stage_summ_AIC <- function(p, R, betaZX, betaZY, se_betaZY, n1, n2, gamma_hat) {
  ZTZ <- R
  ZTY <- diag(ZTZ) * betaZY
  
  YTY_vec <- (n2 - 1) * diag(ZTZ) * se_betaZY^2 + ZTY * betaZY
  YTY <- median(YTY_vec)
  
  calculate_aic <- function(i, ZTY, ZTZ, gamma_hat, n2, YTY) {
    test11 <- diag(0, nrow = p, ncol = p)
    diag(test11)[i] <- 1
    W1 <- cbind(test11, gamma_hat)
    
    solve.W1 <- t(W1) %*% ZTZ %*% W1
    non0 <- rowSums(solve.W1 != 0)
    solve.W1 <- solve.W1[non0 > 0, ]
    non0 <- colSums(solve.W1 != 0)
    solve.W1 <- solve.W1[, non0 > 0]
    
    W1 <- W1[, colSums(W1 != 0) > 0]
    
    beta1 <- ginv(solve.W1) %*% t(W1) %*% ZTY
    n2 * log(YTY - t(beta1) %*% t(W1) %*% ZTY) + 2 * sum(diag(test11))
  }
  
  # Sequential computation of AIC
  testaic <- sapply(1:p, function(i) calculate_aic(i, ZTY, ZTZ, gamma_hat, n2, YTY))
  
  whichIV <- numeric(p)
  whichIV[1] <- which.min(testaic)
  AICtest <- numeric(p)
  AICtest[1] <- testaic[whichIV[1]]
  
  for (j in 2:p) {
    testaic <- sapply(1:p, function(i) {
      test11 <- diag(0, nrow = p, ncol = p)
      diag(test11)[whichIV[1:(j - 1)]] <- 1
      diag(test11)[i] <- 1
      
      W1 <- cbind(test11, gamma_hat)
      solve.W1 <- t(W1) %*% ZTZ %*% W1
      non0 <- rowSums(solve.W1 != 0)
      solve.W1 <- solve.W1[non0 > 0, ]
      non0 <- colSums(solve.W1 != 0)
      solve.W1 <- solve.W1[, non0 > 0]
      
      W1 <- W1[, colSums(W1 != 0) > 0]
      
      beta1 <- ginv(solve.W1) %*% t(W1) %*% ZTY
      n2 * log(YTY - t(beta1) %*% t(W1) %*% ZTY) + 2 * sum(diag(test11))
    })
    
    whichIV[j] <- which.min(testaic)
    AICtest[j] <- testaic[whichIV[j]]
    
    if ((j > 1 && whichIV[j] == whichIV[j - 1]) || j == p) {
      break
    }
  }
  
  first_stage_IV <- unique(whichIV)
  return(first_stage_IV)
}


first_stage_summ_BIC <- function(p, R, betaZX, betaZY, se_betaZY, n1, n2, gamma_hat) {
  # Precompute ZTZ and ZTY
  ZTZ <- R
  ZTY <- diag(ZTZ) * betaZY
  
  # Vectorize YTY calculation
  YTY_vec <- (n2 - 1) * diag(ZTZ) * se_betaZY^2 + ZTY * betaZY
  YTY <- median(YTY_vec)
  
  calculate_bic <- function(i, ZTY, ZTZ, gamma_hat, n2, YTY) {
    test11 <- diag(0, nrow = p, ncol = p)
    diag(test11)[i] <- 1
    W1 <- cbind(test11, gamma_hat)
    
    solve.W1 <- t(W1) %*% ZTZ %*% W1
    non0 <- rowSums(solve.W1 != 0)
    solve.W1 <- solve.W1[non0 > 0, ]
    non0 <- colSums(solve.W1 != 0)
    solve.W1 <- solve.W1[, non0 > 0]
    
    W1 <- W1[, colSums(W1 != 0) > 0]
    
    beta1 <- ginv(solve.W1) %*% t(W1) %*% ZTY
    n2 * log(YTY - t(beta1) %*% t(W1) %*% ZTY) + log(n2) * sum(diag(test11))
  }
  
  # Sequential computation of BIC
  testbic <- sapply(1:p, function(i) calculate_bic(i, ZTY, ZTZ, gamma_hat, n2, YTY))
  
  whichIV <- numeric(p)
  whichIV[1] <- which.min(testbic)
  BICtest <- numeric(p)
  BICtest[1] <- testbic[whichIV[1]]
  
  for (j in 2:p) {
    testbic <- sapply(1:p, function(i) {
      test11 <- diag(0, nrow = p, ncol = p)
      diag(test11)[whichIV[1:(j - 1)]] <- 1
      diag(test11)[i] <- 1
      
      W1 <- cbind(test11, gamma_hat)
      solve.W1 <- t(W1) %*% ZTZ %*% W1
      non0 <- rowSums(solve.W1 != 0)
      solve.W1 <- solve.W1[non0 > 0, ]
      non0 <- colSums(solve.W1 != 0)
      solve.W1 <- solve.W1[, non0 > 0]
      
      W1 <- W1[, colSums(W1 != 0) > 0]
      
      beta1 <- ginv(solve.W1) %*% t(W1) %*% ZTY
      n2 * log(YTY - t(beta1) %*% t(W1) %*% ZTY) + log(n2) * sum(diag(test11))
    })
    
    whichIV[j] <- which.min(testbic)
    BICtest[j] <- testbic[whichIV[j]]
    
    if ((j > 1 && whichIV[j] == whichIV[j - 1]) || j == p) {
      break
    }
  }
  
  first_stage_IV <- unique(whichIV)
  return(first_stage_IV)
}
