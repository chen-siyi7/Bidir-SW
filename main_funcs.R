first_stage_summ_BIC = function(p, R, betaZX, betaZY, se_betaZY, n1, n2, gamma_hat){
  ZTZ = R
  ZTY = matrix(diag(ZTZ) * betaZY, ncol = 1);
  YTY = NULL
  for(SNP in 1:p){
    YTY[SNP] = (n2-1)*ZTZ[SNP, SNP]*(se_betaZY^2)[SNP] + ZTY[SNP]*betaZY[SNP]
  }
  YTY = median(YTY)
  
  testbic = NULL
  for (i in 1:p){
    test11 = diag(x = 0, nrow = p, ncol = p, names = T)
    diag(test11)[i] = 1
    W1 = cbind(test11, gamma_hat)
    solve.W1 = t(W1) %*% ZTZ %*% W1
    non0 = as.numeric(rowSums(solve.W1 != 0))
    solve.W1 = solve.W1[which(non0>0),]
    non0 = as.numeric(colSums(solve.W1 != 0))
    solve.W1 = solve.W1[,which(non0>0)]
    non0 = as.numeric(colSums(W1 != 0))
    W1 = W1[,which(non0>0)]
    beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
    testbic[i] = n2*log(YTY - t(beta1)%*%t(W1)%*%ZTY)+log(n2)*sum(diag(test11))
  }
  
  # one by one 
  whichIV = NULL
  whichIV[1] = which.min(testbic)
  BICtest = NULL
  BICtest[1] = testbic[which.min(testbic)]
  for (j in 2:p){
    testbic = NULL
    for (i in 1:p){
      test11 = diag(x = 0, nrow = p, ncol = p, names = T)
      diag(test11)[whichIV] = 1
      diag(test11)[i] = 1
      W1 = cbind(test11, gamma_hat)
      solve.W1 = t(W1) %*% ZTZ %*% W1
      non0 = as.numeric(rowSums(solve.W1 != 0))
      solve.W1 = solve.W1[which(non0>0),]
      non0 = as.numeric(colSums(solve.W1 != 0))
      solve.W1 = solve.W1[,which(non0>0)]
      non0 = as.numeric(colSums(W1 != 0))
      W1 = W1[,which(non0>0)]
      beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
      testbic[i] = n2*log(YTY - t(beta1)%*%t(W1)%*%ZTY)+log(n2)*sum(diag(test11))
    }
    whichIV[j] = which.min(testbic)
    BICtest[j] = testbic[which.min(testbic)]
    if ((j > 1 && length(whichIV) >= j && !is.na(whichIV[j]) && !is.na(whichIV[j - 1]) && whichIV[j] == whichIV[j - 1]) || j == p) {
      break
    }
  }
  first.stage.IV = unique(whichIV)
  return(first.stage.IV)
}

first_stage_summ_AIC = function(p, R, betaZX, betaZY, se_betaZY, n1, n2, gamma_hat){
  ZTZ = R
  ZTY = matrix(diag(ZTZ) * betaZY, ncol = 1);
  YTY = NULL
  for(SNP in 1:p){
    YTY[SNP] = (n2-1)*ZTZ[SNP, SNP]*(se_betaZY^2)[SNP] + ZTY[SNP]*betaZY[SNP]
  }
  YTY = median(YTY)
  
  testaic = NULL
  for (i in 1:p){
    test11 = diag(x = 0, nrow = p, ncol = p, names = T)
    diag(test11)[i] = 1
    W1 = cbind(test11, gamma_hat)
    solve.W1 = t(W1) %*% ZTZ %*% W1
    non0 = as.numeric(rowSums(solve.W1 != 0))
    solve.W1 = solve.W1[which(non0>0),]
    non0 = as.numeric(colSums(solve.W1 != 0))
    solve.W1 = solve.W1[,which(non0>0)]
    non0 = as.numeric(colSums(W1 != 0))
    W1 = W1[,which(non0>0)]
    beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
    testaic[i] = n2*log(YTY - t(beta1)%*%t(W1)%*%ZTY)+2*sum(diag(test11))
  }
  
  # one by one 
  whichIV = NULL
  whichIV[1] = which.min(testaic)
  AICtest = NULL
  AICtest[1] = testaic[which.min(testaic)]
  for (j in 2:p){
    testaic = NULL
    for (i in 1:p){
      test11 = diag(x = 0, nrow = p, ncol = p, names = T)
      diag(test11)[whichIV] = 1
      diag(test11)[i] = 1
      W1 = cbind(test11, gamma_hat)
      solve.W1 = t(W1) %*% ZTZ %*% W1
      non0 = as.numeric(rowSums(solve.W1 != 0))
      solve.W1 = solve.W1[which(non0>0),]
      non0 = as.numeric(colSums(solve.W1 != 0))
      solve.W1 = solve.W1[,which(non0>0)]
      non0 = as.numeric(colSums(W1 != 0))
      W1 = W1[,which(non0>0)]
      beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
      testaic[i] = n2*log(YTY - t(beta1)%*%t(W1)%*%ZTY)+2*sum(diag(test11))
    }
    whichIV[j] = which.min(testaic)
    AICtest[j] = testaic[which.min(testaic)]
    if ((j > 1 && length(whichIV) >= j && !is.na(whichIV[j]) && !is.na(whichIV[j - 1]) && whichIV[j] == whichIV[j - 1]) || j == p) {
      break
    }
  }
  first.stage.IV = unique(whichIV)
  return(first.stage.IV)
}
