stage2 = function(p, R, betaZX, betaZY, se_betaZY, n1, n2, gamma_hat, sigma2_hat){
  p = p
  ZTZ = R
  ZTY = matrix(diag(ZTZ) * betaZY, ncol = 1);
  YTY = NULL
  for(SNP in 1:p){
    YTY[SNP] = (n2-1)*ZTZ[SNP, SNP]*(se_betaZY^2)[SNP] + ZTY[SNP]*betaZY[SNP]
  }
  YTY = median(YTY)
  # result
  W1 = gamma_hat
  solve.W1 = t(W1) %*% ZTZ %*% W1
  beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
  beta_est = as.numeric(tail(beta1, n = 1))
  C = ginv(solve.W1)%*%(t(W1)%*%ZTZ)
  sigma_u2 = as.numeric((YTY - t(beta1)%*%t(W1)%*%ZTY))
  Varbeta = ginv(solve.W1*n2)*sigma_u2+beta_est^2*sigma2_hat*ginv(solve.W1*n2)*n2/n1
  Varbeta = tail(diag(Varbeta), n = 1)
  beta_se = as.numeric(sqrt(Varbeta))
  Z = beta_est/beta_se
  Pvalue = 2*pnorm(-abs(Z), 0, 1)
  my_list = list('Estimate' = beta_est, 'SE' = beta_se, 'Pvalue' = Pvalue)
  return(my_list)
}

stage2_egger = function(p, R, betaZX, betaZY, se_betaZY, n1, n2, gamma_hat, sigma2_hat){
  p = p
  ZTZ = R
  ZTY = matrix(diag(ZTZ) * betaZY, ncol = 1);
  YTY = NULL
  for(SNP in 1:p){
    YTY[SNP] = (n2-1)*ZTZ[SNP, SNP]*(se_betaZY^2)[SNP] + ZTY[SNP]*betaZY[SNP]
  }
  YTY = median(YTY)
  # result
  test11 = diag(x = 1, nrow = p, ncol = p, names = T)
  W1 = cbind(test11, gamma_hat)
  solve.W1 = t(W1) %*% ZTZ %*% W1
  beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
  beta_est = as.numeric(tail(beta1, n = 1))
  C = ginv(solve.W1)%*%(t(W1)%*%ZTZ)
  sigma_u2 = as.numeric((YTY - t(beta1)%*%t(W1)%*%ZTY))
  Varbeta = ginv(solve.W1*n2)*sigma_u2+beta_est^2*sigma2_hat*ginv(solve.W1*n2)*n2/n1
  Varbeta = tail(diag(Varbeta), n = 1)
  beta_se = as.numeric(sqrt(Varbeta))
  Z = beta_est/beta_se
  Pvalue = 2*pnorm(-abs(Z), 0, 1)
  my_list = list('Estimate' = beta_est, 'SE' = beta_se, 'Pvalue' = Pvalue)
  return(my_list)
}

