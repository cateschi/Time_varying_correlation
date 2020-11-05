## Packages needed to run the script ##

library(ucminf)
library(magic)      # to create block-diagonal matrices
library(numDeriv)      # to compute the Hessian and the gradient numerically
library(ggplot2)      # for plots
library(cowplot)      # for plots
library(matrixcalc)      # for vec operator
library(MASS)      # to draw from a multivariate normal distribution
library(matrixStats)      # to get quantiles out of matrix
library(Rcpp)
library(RcppArmadillo)
library(zoo)
library(PerformanceAnalytics)
library(optimr)



## Functions ##

# Function that creates the cubic splines weights
Weights <- function(len, knots){
  
  n <- len   
  time <- c(1:n)      
  k <- length(knots)     
  h <- diff(knots)      
  
  lambda <- rep(NA,(k-2))
  for (j in 1:(k-2)){
    lambda [j] <- h[j+1]/(h[j]+h[j+1])
  }
  
  pi.zero <- 0
  pi.k <- 0
  
  Lambda <- matrix(NA, k, k)     
  for (j in 1:k){
    if (j == 1){
      Lambda[j,] <- c(c(2, -2*pi.zero), rep(0,(k-2)))
    } else if (j > 1 && j < k){
      Lambda[j,] <- c(rep(0,(j-2)), c((1-lambda[j-1]), 2, lambda[j-1]), rep(0,(k-length(rep(0,(j-2)))-length(c((1-lambda[j-1]), 2, lambda[j-1])))))
    } else if (j == k){
      Lambda[j,] <- c(rep(0,(k-2)), c(-2*pi.zero, 2))
    }
  }
  
  Theta <- matrix(NA, k, k)     
  for (j in 1:k){
    if (j == 1){
      Theta[j,] <- rep(0,k)
    } else if (j > 1 && j < k){
      Theta[j,] <- c(rep(0,(j-2)), c(6/(h[j-1]*(h[j-1]+h[j])), -6/(h[j-1]*h[j]), 6/(h[j]*(h[j-1]+h[j]))), rep(0,(k-length(rep(0,(j-2)))-length(c((1-lambda[j-1]), 2, lambda[j-1])))))
    } else if (j == k){
      Theta[j,] <- rep(0,k)
    }
  }
  
  P <- matrix(0, n, k)    
  for (i in 1:n){      
    for (j in 2:k){     
      if (i >= (knots[j-1]) && i <= knots[j]){
        P[i,(j-1)] <- (knots[j] - i)/(6*h[j-1])*((knots[j] - i)^2 - h[j-1]^2)
        P[i,j] <- (i - knots[j-1])/(6*h[j-1])*((i - knots[j-1])^2 - h[j-1]^2)
      }
    }
  }
  
  Q <- matrix(0, n, k)     
  for (i in 1:n){      
    for (j in 2:k){     
      if (i >= (knots[j-1]) && i <= knots[j]){
        Q[i,(j-1)] <- (knots[j] - i)/h[j-1]
        Q[i,j] <- (i - knots[j-1])/h[j-1]
      }
    }
  }
  
  W <- P%*%solve(Lambda)%*%Theta + Q      
  
  return(list(W=W, P=P, Q=Q, Lambda=Lambda, Theta=Theta))
}


link <- function(x){
  y <- x/sqrt(1+x^2)
  return(y)
}


KF_CC_known_corr <- function(par,y,se,opti,outofsample,parP10,nstates,hyper_tan,d,gamma_draw){
  len <- length(y[1,])
  sigma_Ry <- par[1]
  sigma_omegay <- par[2]
  sigma_lambda <- par[3]
  sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
  sigma_Rx <- par[9]
  sigma_omegax <- par[10]
  x10 <- rep(0,nstates)
  Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
  Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
  delta <- par[12]
  P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
  Pttm1[[1]] <- P10
  xtt <- matrix(0,nstates,(len))
  st_for <- matrix(0,nrow(y),len)
  xttm1 <- matrix(0,nstates,(len+1))
  xttm1[,1] <- x10
  H <- adiag(diag(0,5,5), exp(2*par[11]))
  
  #Bulid T:
  Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
  C <- array(0,dim=c(2,2,5))
  for (l in 1:5){
    C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
  }
  Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
  ncol(Tyomega)
  nrow(Tyomega)
  Tylambda <- diag(4)
  TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
  delta <- delta
  TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
  TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
  Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
  Tx <- adiag(Tymu, Tyomega)
  Tmatrix <- adiag(Ty, Tx)
  
  R <- diag(1,nstates,nstates)
  D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11))
  
  #initialization of loglikelihood
  logl <- 0
  
  #Start of KF recursions
  for (i in 1:len){ 
    
    #Bulid Z:
    Zy <- c(1,0)
    Zy <- rep(Zy,6)
    Zy <- c(Zy,1)
    Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
    Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
    Zy <- cbind(Zy, diag(as.numeric(se[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
    Zx <- c(1,0)
    Zx <- rep(Zx,6)
    Zx <- c(Zx,1)
    Zx <- rbind(Zx)
    Z <- adiag(Zy,Zx)
    ncol(Z)
    nrow(Z)
    
    epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
    Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + H
    #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
    if ((NaN %in% Fmatrix)==T){
      logl<- -P10[1]
    } else {
      svdFmatrix <- svd(Fmatrix)
      Fmatrix_inv <- svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)
      if (opti==F){
        Bmatrix <- chol(Fmatrix_inv)
        st_for[,i] <- Bmatrix%*%epshatoutofsample
      }
      Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv #kalman gain
      xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%epshatoutofsample #compute x_{t|t}
      epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
      Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%Z%*%Pttm1[[i]] #compute P_{t|t}
      
      if (hyper_tan==T){
        R[32,2] <- tanh(gamma_draw[i]) 
        R[2,32] <- tanh(gamma_draw[i]) 
       } else {
        R[32,2] <- link(gamma_draw[i])      
        R[2,32] <- link(gamma_draw[i]) 
      }
      Q <- D%*%R%*%D 
      
      Pttm1[[i+1]] <- Tmatrix%*%Pttm1[[i]]%*%t(Tmatrix-Kg%*%Z)+Q 
      xttm1[,i+1] <- Tmatrix%*%xttm1[,i] + Kg%*%epshatoutofsample
      
      #The optimization criterion
      if (outofsample) {
        if (i <= d ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > d ){ # diffuse log likelihood
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix_inv%*%epshatoutofsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      } else {
        if (i <= d ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > d ){ # diffuse log likelihood
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%Fmatrix_inv%*%epshatinsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      }
    }
  }
  if (opti) {
    return(-logl)
  }
  else {
    return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt,st_for=st_for))
  }
}


KF_CC_splines <- function(par,y,se,opti,outofsample,parP10,nstates,k,W,restricted,hyper_tan,d){
  len <- length(y[1,])
  sigma_Ry <- par[1]
  sigma_omegay <- par[2]
  sigma_lambda <- par[3]
  sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
  sigma_Rx <- par[9]
  sigma_omegax <- par[10]
  x10 <- rep(0,nstates)
  Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
  Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
  delta <- par[12]
  P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
  Pttm1[[1]] <- P10
  xtt <- matrix(0,nstates,(len))
  st_for <- matrix(0,nrow(y),len)
  xttm1 <- matrix(0,nstates,(len+1))
  xttm1[,1] <- x10
  H <- adiag(diag(0,5,5), exp(2*par[11]))
  
  #Bulid T:
  Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
  C <- array(0,dim=c(2,2,5))
  for (l in 1:5){
    C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
  }
  Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
  ncol(Tyomega)
  nrow(Tyomega)
  Tylambda <- diag(4)
  TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
  delta <- delta
  TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
  TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
  Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
  Tx <- adiag(Tymu, Tyomega)
  Tmatrix <- adiag(Ty, Tx)
  
  R <- diag(1,nstates,nstates)
  D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11))
  
  #initialization of loglikelihood
  logl <- 0
  
  #Start of KF recursions
  for (i in 1:len){ 
    
    #Bulid Z:
    Zy <- c(1,0)
    Zy <- rep(Zy,6)
    Zy <- c(Zy,1)
    Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
    Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
    Zy <- cbind(Zy, diag(as.numeric(se[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
    Zx <- c(1,0)
    Zx <- rep(Zx,6)
    Zx <- c(Zx,1)
    Zx <- rbind(Zx)
    Z <- adiag(Zy,Zx)
    ncol(Z)
    nrow(Z)
    
    epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
    Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + H
    #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
    if ((NaN %in% Fmatrix)==T){
      logl<- -P10[1]
    } else {
      svdFmatrix <- svd(Fmatrix)
      Fmatrix_inv <- svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)
      if (opti==F){
        Bmatrix <- chol(Fmatrix_inv)
        st_for[,i] <- Bmatrix%*%epshatoutofsample
      }
      Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv #kalman gain
      xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%epshatoutofsample #compute x_{t|t}
      epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
      Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%Z%*%Pttm1[[i]] #compute P_{t|t}
      
      if (hyper_tan==T){  
        if (restricted==T){
          R[32,2] <- tanh(par[length(par)-k+1])      # time constant correlation
          R[2,32] <- tanh(par[length(par)-k+1]) 
        } else {
          R[32,2] <- tanh(W[i,1:k]%*%par[(length(par)-k+1):length(par)])      # time varying correlation
          R[2,32] <- tanh(W[i,1:k]%*%par[(length(par)-k+1):length(par)]) 
        }
      } else {
        if (restricted==T){
          R[32,2] <- link(par[length(par)-k+1])      # time constant correlation
          R[2,32] <- link(par[length(par)-k+1]) 
        } else {
          R[32,2] <- link(W[i,1:k]%*%par[(length(par)-k+1):length(par)])      # time varying correlation
          R[2,32] <- link(W[i,1:k]%*%par[(length(par)-k+1):length(par)]) 
        }
      }
      Q <- D%*%R%*%D 
      
      Pttm1[[i+1]] <- Tmatrix%*%Pttm1[[i]]%*%t(Tmatrix-Kg%*%Z)+Q 
      xttm1[,i+1] <- Tmatrix%*%xttm1[,i] + Kg%*%epshatoutofsample
      
      #The optimization criterion
      if (outofsample) {
        if (i <= d ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > d ){ # diffuse log likelihood
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix_inv%*%epshatoutofsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      } else {
        if (i <= d ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > d ){ # diffuse log likelihood
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%Fmatrix_inv%*%epshatinsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      }
    }
  }
  if (opti) {
    return(-logl)
  }
  else {
    return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt,st_for=st_for))
  }
}


KF_CC_const <- function(par,y,se,opti,outofsample,parP10,nstates,hyper_tan){
  len <- length(y[1,])
  sigma_Ry <- par[1]
  sigma_omegay <- par[2]
  sigma_lambda <- par[3]
  sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
  sigma_Rx <- par[9]
  sigma_omegax <- par[10]
  x10 <- rep(0,nstates)
  Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
  Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
  delta <- par[13]
  P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
  Pttm1[[1]] <- P10
  xtt <- matrix(0,nstates,(len))
  st_for <- matrix(0,nrow(y),len)
  xttm1 <- matrix(0,nstates,(len+1))
  xttm1[,1] <- x10
  R <- diag(1,nstates,nstates)
  D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11))
  if (hyper_tan==T){
    R[32,2] <- tanh(par[11])
    R[2,32] <- tanh(par[11])
  } else {
    R[32,2] <- link(par[11])
    R[2,32] <- link(par[11])
  }
  Q <- D%*%R%*%D
  H <- adiag(diag(0,5,5), exp(2*par[12]))
  
  #Bulid T:
  Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
  C <- array(0,dim=c(2,2,5))
  for (l in 1:5){
    C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
  }
  Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
  ncol(Tyomega)
  nrow(Tyomega)
  Tylambda <- diag(4)
  TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
  delta <- delta
  TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
  TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
  Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
  Tx <- adiag(Tymu, Tyomega)
  Tmatrix <- adiag(Ty, Tx)
  
  #initialization of loglikelihood
  logl <- 0
  
  #Start of KF recursions
  for (i in 1:len){ 
    
    #Bulid Z:
    Zy <- c(1,0)
    Zy <- rep(Zy,6)
    Zy <- c(Zy,1)
    Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
    Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
    Zy <- cbind(Zy, diag(as.numeric(se[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
    Zx <- c(1,0)
    Zx <- rep(Zx,6)
    Zx <- c(Zx,1)
    Zx <- rbind(Zx)
    Z <- adiag(Zy,Zx)
    ncol(Z)
    nrow(Z)
    
    epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
    Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + H
    #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
    if ((NaN %in% Fmatrix)==T){
      logl<- -P10[1]
    } else {
      svdFmatrix <- svd(Fmatrix)
      Fmatrix_inv <- svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)
      if (opti==F){
        Bmatrix <- chol(Fmatrix_inv)
        st_for[,i] <- Bmatrix%*%epshatoutofsample
      }
      Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv #kalman gain
      xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%epshatoutofsample #compute x_{t|t}
      epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
      Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%Z%*%Pttm1[[i]] #compute P_{t|t}
      Pttm1[[i+1]] <- Tmatrix%*%Pttm1[[i]]%*%t(Tmatrix-Kg%*%Z)+Q 
      xttm1[,i+1] <- Tmatrix%*%xttm1[,i] + Kg%*%epshatoutofsample
      
      #The optimization criterion
      if (outofsample) {
        if (i <= (30-13) ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > (30-13) ){ # diffuse log likelihood
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix_inv%*%epshatoutofsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      } else {
        if (i <= (30-13) ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > (30-13) ){ # diffuse log likelihood
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%Fmatrix_inv%*%epshatinsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      }
    }
  }
  if (opti) {
    return(-logl)
  }
  else {
    return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt,st_for=st_for))
  }
}


IndInf_CC <- function(par, sim_h, beta_hat, len, eps_sim, VARgamma, opti,outofsample,parP10,nstates,d,k,W,restricted,hyper_tan,se,states_noerr,nvar,init_val_CC){
  
  sigma_Ry <- par[1]
  sigma_omegay <- par[2]
  sigma_lambda <- par[3]
  sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
  sigma_Rx <- par[9]
  sigma_omegax <- par[10]
  
  betatilde <- matrix(NA, nrow=sim_h, ncol=length(par))
  
  possibleError <- tryCatch(for (h in 1:sim_h){
    
    # generate TV corr according to a random walk model:
    gammaII <- cumsum(exp(par[length(par)])*eps_sim[2,(h*len-len+1):(h*len)])
    if (hyper_tan == T){
      rhoII <- tanh(gammaII)
    } else {
      rhoII <- link(gammaII)
    }
    
    # generate Innovations of states:
    Q <- lapply(seq_len(len), function(X) matrix(0,nstates-length(states_noerr),nstates-length(states_noerr)))
    R <- lapply(seq_len(len), function(X) diag(1,nstates-length(states_noerr),nstates-length(states_noerr)))
    D <- adiag(exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, exp(sigma_Rx), exp(sigma_omegax)*diag(11))
    M <- lapply(seq_len(len), function(X) matrix(0,nstates-length(states_noerr),nstates-length(states_noerr)))
    eta <- eps_sim[3:nrow(eps_sim),(h*len-len+1):(h*len)]      # innovations of the state equation
    for (i in 1:len){
      R[[i]][1,22] <- rhoII[i]
      R[[i]][22,1] <- rhoII[i]
      Q[[i]] <- D%*%R[[i]]%*%D
      M[[i]] <- t(chol(D)%*%chol(R[[i]])%*%chol(D))
      eta[,i] <- M[[i]]%*%eta[,i]
    }
    
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- par[12]
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(Tymu, Tyomega)
    Tmatrix <- adiag(Ty, Tx)
    
    
    # generate state variables
    Rsel <- diag(1,nstates)      # selection matrix of the innovations in the transition equation
    Rsel <- Rsel[,-states_noerr]
    alpha <- matrix(0,nstates,len)      # state vector
    alpha[,1] <- Rsel%*%eta[,1]
    for (i in 2:len){
      alpha[,i] <- Tmatrix%*%alpha[,i-1] + Rsel%*%eta[,i]
    }
    #for (j in 1:nrow(alpha)){plot(alpha[j,], type="l")}
    
    
    # Generate the innovations in the observation equation
    nu <- exp(par[11])*eps_sim[1,(h*len-len+1):(h*len)]
    
    Msel <- diag(1,nvar)      # selection matrix of the innovations in the observation equation
    Msel <- Msel[,-c(1:5)]
    
    # Generate the observable variables - the final model 
    y <- matrix(NA,nvar,len)
    for (i in 1:len){
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(se[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- c(1,0)
      Zx <- rep(Zx,6)
      Zx <- c(Zx,1)
      Zx <- rbind(Zx)
      Z <- adiag(Zy,Zx)
      
      y[,i] <- Z%*%alpha[,i] + Msel%*%as.matrix(nu[i])
    }
    #for (j in 1:nrow(y)){plot(y[j,], type="l")}
    
    
    init_valII <- c(beta_hat,rep(0,length(knots)))
    
    objoptII <- ucminf_rcpp_splines(init_valII,y=y,se=as.matrix(se),opti=opti,outofsample=outofsample,
                                    parP10=parP10,nstates=nstates, d=d,k=k, W=W, hyper_tan=hyper_tan, restricted=restricted, 
                                    control=list(gradstep = c(1e-2, 1e-3)))
    
    betatilde[h,] <- c(objoptII$par[-c((length(objoptII$par)-length(knots)+1):length(objoptII$par))], var(objoptII$par[13:length(objoptII$par)]))
  } , error = function(e) e)
  if(inherits(possibleError, "error")) {
    opt.fun <- 1000000000
  } else {
    opt.fun <- t(c(beta_hat, VARgamma) - colMeans(betatilde))%*%(c(beta_hat, VARgamma) - colMeans(betatilde))
  }
  
  
  #if (abs(theta.opt) > 1){opt.fun <- 100000}      # push away parameter if explosive MA process
  
  return(opt.fun)
}




#### Set working directory ####

si <- Sys.info()
if (si[7] == "Caterina") {
  base_folder <- "C:/Users/Caterina/Dropbox/"
} else if (si[7] == "c.schiavoni") {
  base_folder <- "C:/Users/C.Schiavoni/Dropbox/"
} else if (si[7] == "C.Schiavoni") {
  base_folder <- "C:/Users/C.Schiavoni/Dropbox/"
}
folder <- paste(base_folder,"Time_varying_correlation/Empirics/",sep="")
setwd(folder)



#### Upload cpp function

Cpp_IndInf_ssm_levels <- sourceCpp("estimation_CC_IndInf_BF.cpp")



#### Upload data ####

sep <- FALSE     # if TRUE, upload data until Sep 2020, otherwise until March 2020

if (sep == T){
  
  LFS <- read_excel("data/InputseriesCaterina_Sep2020.xlsx", col_names = FALSE)
  
  CC <- read_excel("data/ww_uitkeringen_Sep2020.xls", col_names = FALSE)
  
  se <- LFS[-c(1:37),-c(1:9)]
  
} else {

  LFS <- read_excel("data/InputseriesCaterina_Mar2020.xlsx", col_names = FALSE)
  
  CC <- read_excel("data/ww_uitkeringen_Mar2020.xls", col_names = FALSE)
  
  se <- read_excel("data/SFInputSeries_Mar2020.xlsx", col_names = FALSE)
  se <- as.data.frame(se[-c(1:37),])

}

LFS <- LFS[-c(1:37),]
waves <- LFS[,c(3:7)]
waves <- as.data.frame(waves)
waves <- lapply(waves, as.numeric)
waves <- as.data.frame(waves)
waves <- as.matrix(waves)

len_m <- nrow(waves)

CC <- CC[-c(1:36),1]
CC <- as.data.frame(CC)
CC <- lapply(CC, as.numeric)
CC <- as.data.frame(CC)
CC <- as.matrix(CC)

se <- lapply(se, as.numeric)
se <- as.data.frame(se)
se <- as.matrix(se)




#### LFS and CC ####

len_m <- nrow(waves)

y <- matrix(0,6,len_m)
for (j in 1:5){
  y[j,] <- waves[1:len_m,j]
}
y[6,] <- CC[1:len_m,1]
#ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,],y[6,]),col=c("black","black","black","black","black","red"))

hyper_tan <- T      # if TRUE, use the hyperbolic tangent as link function for the TV correlation

init_val_CC <- c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                 log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),
                 log(3000),log(0.02),0,log(1000), 0.21)



#### Time constant correlation ####

objopt_CC_const <- ucminf(par=init_val_CC,
                           KF_CC_const,y=y,se=se,opti=T,outofsample=T,parP10=1000000000000,nstates=43, hyper_tan=hyper_tan,  hessian=2,
                           control=list(gradstep = c(1e-2, 1e-3)))

par_biv_const <- objopt_CC_const$par
#exp(par_biv_const)/(1-par_biv_const[13]^2) #to get the estimates for the standard deviation of the sampling errors

obj <- KF_CC_const(par_biv_const,y=y,se=se,opti=F,outofsample=T,parP10=1000000000000,nstates=43, hyper_tan=hyper_tan)

P_L_biv <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
P_R_biv <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
P_theta_biv <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]

L_biv <- obj$xtt[1,(31-13):len_m]
R_biv <- obj$xtt[2,(31-13):len_m]
season_y_est_tt <- obj$xtt[3,]+obj$xtt[5,]+obj$xtt[7,]+obj$xtt[9,]+obj$xtt[11,]+obj$xtt[13,]      #estimated seasonal components
theta_biv <- obj$xtt[1,] + season_y_est_tt 
theta_biv <- theta_biv[(31-13):len_m]

logl_const <- -obj$logl

res_const <- obj$st_for      # standardised residuals


#plot(c(corr_2010,corr_2011,corr_2012,corr_2013,corr_2014,corr_2020), ylab=" ", xlab=" ",xaxt = 'n', pch=20, cex=2.5)
#axis(1, at=c(1:6), labels=c("2010", "2011", "2012", "2013", "2014", "2020"))



#### Cubic splines ####


knotFeb <- F      # if TRUE knot at February 2015

AIC <- rep(NA,2)
BIC <- rep(NA,2)

for (choice_k in 1:2){
  
  if (choice_k==1){ 
    knots <- as.numeric(quantile(c(1:len_m), probs = seq(0, 1, 0.25)))
  } else {
    knots <- as.numeric(quantile(c(1:len_m), probs = seq(0, 1, 1/7)))
  }
  
  W <- Weights(len=len_m,knots=knots)$W
  
  
  # Kalman filter estimation
  
  # objopt_CC_splines <- ucminf(par=c(init_val_CC[-11],rep(0,length(knots))),
  #                             KF_slopes_CC_splines,y=y,se=se,opti=T,outofsample=T,parP10=1000000000000,nstates=43, 
  #                             k=length(knots), restricted=F, hyper_tan=hyper_tan, W=W, d=30-13, hessian=2 , control=list(grad="central", gradstep = c(1e-2, 1e-3), trace=T))
  
  objopt_CC_splines <- ucminf_rcpp_splines(init_val=c(init_val_CC[-11],rep(0,length(knots))),y=y,se=as.matrix(se),opti=T,outofsample=T,
                                           parP10=1000000000000,nstates=43, d=30-13,k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, 
                                           control=list(gradstep = c(1e-2, 1e-3)))
  
  
  AIC[choice_k] <- 1/len_m*(2*objopt_CC_splines$value + 2*(30-13 + length(c(init_val_CC[-11],rep(0,length(knots))))))
  BIC[choice_k] <- 1/len_m*(2*objopt_CC_splines$value + log(len_m)*(30-13 + length(c(init_val_CC[-11],rep(0,length(knots))))))
  
}

if (which.min(BIC) == 1){
  knots <- as.numeric(quantile(c(1:len_m), probs = seq(0, 1, 0.25)))
  
  W <- Weights(len=len_m,knots=knots)$W
  
  start_time_CC_splines <- proc.time()      # Start the clock
  objopt_CC_splines <- ucminf_rcpp_splines(init_val=c(init_val_CC[-11],rep(0,length(knots))),y=y,se=as.matrix(se),opti=T,outofsample=T,
                                           parP10=1000000000000,nstates=43, d=30-13,k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, 
                                           control=list(gradstep = c(1e-2, 1e-3)))
  el_time_CC_splines <- proc.time() - start_time_CC_splines      # Stop the clock
}


if (knotFeb==T){
  knots <- as.numeric(c(quantile(c(1:134), probs = seq(0, 1, 0.4)), 134,len_m))
  
  W <- Weights(len=len_m,knots=knots)$W
  
  objopt_CC_splines <- ucminf_rcpp_splines(init_val=c(init_val_CC[-11],rep(0,length(knots))),y=y,se=as.matrix(se),opti=T,outofsample=T,
                                           parP10=1000000000000,nstates=43, d=30-13,k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, 
                                           control=list(gradstep = c(1e-2, 1e-3)))
}


#saveRDS(objopt_CC_splines, file.path("insample_objopt_CC_splines.csv"))

#objopt_CC_splines <- readRDS(file.path("insample_objopt_CC_splines.csv"))

par_biv_splines <- objopt_CC_splines$par
#exp(par_biv_splines)/(1-par_biv_splines[12]^2) #to get the estimates for the standard deviation of the sampling errors

obj <- KF_CC_splines(objopt_CC_splines$par,y=y,se=se,opti=F,outofsample=T,parP10=1000000000000,nstates=43,
                            k=length(knots), restricted=F, W=W, d=30-13, hyper_tan=hyper_tan)

P_L_biv_splines <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
P_R_biv_splines <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
P_theta_biv_splines <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]

L_biv_splines <- obj$xtt[1,(31-13):len_m]
R_biv_splines <- obj$xtt[2,(31-13):len_m]
season_y_est_tt <- obj$xtt[3,]+obj$xtt[5,]+obj$xtt[7,]+obj$xtt[9,]+obj$xtt[11,]+obj$xtt[13,]      #estimated seasonal components
theta_biv_splines <- obj$xtt[1,] + season_y_est_tt 
theta_biv_splines <- theta_biv_splines[(31-13):len_m]

logl_splines <- -obj$logl

gamma_hat_splines <- W[,1:length(knots)]%*%par_biv_splines[(length(par_biv_splines)-length(knots)+1):length(par_biv_splines)]
if (hyper_tan==T){
  rho_hat_splines <- tanh(gamma_hat_splines)
} else {
  rho_hat_splines <- link(gamma_hat_splines)
}
plot(rho_hat_splines, type="l")

res_splines <- obj$st_for      # standardised residuals




#### Indirect inference ####

# state variables that do not have an error term in the transition equation
states_noerr <- c(1,23:31)
Rsel <- diag(1,43)      # selection matrix of the innovations in the transition equation
Rsel <- Rsel[,-states_noerr]
Msel <- diag(1,6)      # selection matrix of the innovations in the observation equation
Msel <- Msel[,-c(1:5)]

sim_h <- 5
set.seed(2609)
long.series <- F      # if TRUE, run one cubic splines on TS observations
eps_sim <- t(mvrnorm(n = len_m*sim_h, rep(0,43-length(states_noerr)+2), diag(1,43-length(states_noerr)+2)))  


# find initial value
grid_sigma_gamma <- seq(0.01, 0.15, 0.01)
obj_grid <- rep(NA, length(grid_sigma_gamma))

for (j in 1:length(grid_sigma_gamma)){
  obj_grid[j] <- IndInf_CC(par=c(par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))],log(grid_sigma_gamma[j])), 
                        sim_h=sim_h, beta_hat=par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))], len=len_m, eps_sim=eps_sim, VARgamma = var(par_biv_splines[13:length(par_biv_splines)]),
                        opti=T,outofsample=T,parP10=1000000000000,nstates=43, d=30-13, k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, se=se, states_noerr=states_noerr,nvar=6,init_val_CC=init_val_CC)
}
plot(obj_grid, type="l")
grid_sigma_gamma[which.min(obj_grid)]


IndInf_result <- optimr(par=c(par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))],log(grid_sigma_gamma[which.min(obj_grid)])), IndInf_CC, 
                        sim_h=sim_h, beta_hat=par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))], len=len_m, eps_sim=eps_sim, VARgamma = var(par_biv_splines[13:length(par_biv_splines)]),
                        opti=T,outofsample=T,parP10=1000000000000,nstates=43, d=30-13, k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, se=se, states_noerr=states_noerr,nvar=6,init_val_CC=init_val_CC,
                        control=list(trace=T), method="CG")

#IndInf_result <- readRDS("IndInf_result_CC_S5.csv")
#exp(IndInf_result$par)/(1-IndInf_result$par[12]^2) #to get the estimates for the standard deviation of the sampling errors

# IndInf_result_cpp <- ucminf_rcpp_IndInf(c(par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))],log(0.1)), 
#                         sim_h=sim_h, beta_hat=par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))], eps_sim=eps_sim, VARgamma = var(gamma_hat_splines),
#                         opti=T,outofsample=T,parP10=10000000,d=30-13, W=W, hyper_tan=hyper_tan, restricted=F, se=as.matrix(se), states_noerr=states_noerr,
#                         init_val_CC=init_val_CC, Rsel=Rsel, Msel=Msel, init_valII=c(init_val_CC[-11],rep(0,length(knots))), control=list(gradstep = c(1e-2, 1e-3), trace=T))




#### Bootstrap filter ####

draw_m <- 5000

if (knotFeb==T){
  set.seed(2811)    
} else {
  set.seed(1006)
}
results_BF <- boot_filter_CC_rcpp(draw_m, IndInf_result$par, y, as.matrix(se), nstates=43, hyper_tan, Rsel, states_noerr, init_gamma = par_biv_const[11])

KF_final <- KF_CC_known_corr(par=IndInf_result$par[1:12],y=y[,-len_m],se,opti=F,outofsample=T,parP10=1000000000000,nstates=43,d=30-13,
                             gamma_draw = results_BF$att_BF[2:len_m], hyper_tan = hyper_tan)
att_KF <- KF_final$xtt
Ptt_KF <- KF_final$Ptt

P_L_biv_BF <- unlist(lapply(KF_final$Ptt, function(x) x[1,1]))[(31-13):length(unlist(lapply(KF_final$Ptt, function(x) x[1,1])))]
P_R_biv_BF <- unlist(lapply(KF_final$Ptt, function(x) x[2,2]))[(31-13):length(unlist(lapply(KF_final$Ptt, function(x) x[2,2])))]
P_theta_biv_BF <- unlist(lapply(KF_final$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(31-13):length(unlist(lapply(KF_final$Ptt, function(x) x[2,2])))]

L_biv_BF <- KF_final$xtt[1,(31-13):(len_m-1)]
R_biv_BF <- KF_final$xtt[2,(31-13):(len_m-1)]
season_y_est_tt <- KF_final$xtt[3,]+KF_final$xtt[5,]+KF_final$xtt[7,]+KF_final$xtt[9,]+KF_final$xtt[11,]+KF_final$xtt[13,]      #estimated seasonal components
theta_biv_BF <- KF_final$xtt[1,] + season_y_est_tt 
theta_biv_BF <- theta_biv_BF[(31-13):(len_m-1)]

logl_BF <- -KF_final$logl

res_BF <- KF_final$st_for      # standardised residuals


cv <- 1.96      # critical value for confidence intervals
ts.plot(cbind(att_KF[2,-c(1:(30-13))], att_KF[2,-c(1:(30-13))]+cv*sqrt(unlist(lapply(Ptt_KF, function(x) x[2,2]))[-c(1:(30-13))]), att_KF[2,-c(1:(30-13))]-cv*sqrt(unlist(lapply(Ptt_KF, function(x) x[2,2])))[-c(1:(30-13))]), col=c(rep("black",3)), lty=c(1,2,2))
ts.plot(cbind(att_KF[32,-c(1:(30-13))], att_KF[32,-c(1:(30-13))]+cv*sqrt(unlist(lapply(Ptt_KF, function(x) x[32,32]))[-c(1:(30-13))]), att_KF[32,-c(1:(30-13))]-cv*sqrt(unlist(lapply(Ptt_KF, function(x) x[32,32])))[-c(1:(30-13))]), col=c(rep("black",3)), lty=c(1,2,2))
if (hyper_tan==T){
  ts.plot(cbind(tanh(results_BF$att_BF[-c(1:(30-13),len_m)]), tanh(results_BF$att_BF[-c(1:(30-13),len_m)]+cv*sqrt(results_BF$Ptt_BF[-c(1:(30-13),len_m)])),  tanh(results_BF$att_BF[-c(1:(30-13),len_m)]-cv*sqrt(results_BF$Ptt_BF[-c(1:(30-13),len_m)]))), col=c(rep("black",3)), lty=c(1,2,2))
  ts.plot(cbind(tanh(results_BF$att_BF[-c(1:(30-13),len_m)])), col=c("black"))
} else {
  ts.plot(cbind(link(results_BF$att_BF[-c(1:(30-13),len_m)]), link(results_BF$att_BF[-c(1:(30-13),len_m)]+cv*sqrt(results_BF$Ptt_BF[-c(1:(30-13),len_m)])),  link(results_BF$att_BF-cv*sqrt(results_BF$Ptt_BF[-c(1:(30-13),len_m)]))), col=c(rep("black",3)), lty=c(1,2,2))
}


if (hyper_tan == T){
  rho_hat_BF = tanh(results_BF$att_BF)
} else {
  rho_hat_BF = link(results_BF$att_BF)
}


ts.plot(cbind(obj$xtt[2,-c(1:(30-13),len_m)], att_KF[2,-c(1:(30-13),len_m)]), col=c("black", "red"), lty=c(1,1,2,2))
ts.plot(cbind(obj$xtt[32,-c(1:(30-13),len_m)], att_KF[2,-c(1:(30-13),len_m)]), col=c("black", "red"), lty=c(1,1,2,2))
ts.plot(cbind(rho_hat_splines, rho_hat_BF), col=c("black", "red"), lty=c(1,1,2,2))

plot(results_BF$ESS, type="l")
plot(results_BF$CV[-c(1:4)], type="l")




#### TS plots ####

date <- seq(as.Date("2004/1/1"), by = "month", length.out = len_m-1)

rho_CS_ts <- zoo(x = rho_hat_splines)
rho_BF_ts <- zoo(x = rho_hat_BF)
rho_const_ts <- zoo(x = rep(tanh(par_biv_const[11]), len_m))
gamma_CS_ts <- zoo(x = gamma_hat_splines)
gamma_BF_ts <- zoo(x = results_BF$att_BF)
gamma_const_ts <- zoo(x = rep(par_biv_const[11], len_m))
LFS_ts <- zoo(waves)
CC_ts <- zoo(CC)


fmt <- "%Y"
ix <- seq(1, length(date), 12)
labs <- format(date, fmt)


# plot TV correlation rho

#plot(date, rho_CS_ts[-len_m], lwd=2, main=expression(hat(rho)[t]), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red", ylim=c(-0.5,1))
plot(date, rho_CS_ts[-len_m], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red", ylim=c(-0.5,1))
axis(side = 1, at = date[ix], labels = labs[ix], tcl = -0.7, cex.axis = 0.7)
u <- par("usr")
#plot green rect - notice that x-coordinates are defined by date points
rect(date[49], u[3], date[64], u[4], border = 0, col = "lightgray")      # global financial crisis
#rect(date[92], u[3], date[110], u[4], border = 0, col = "lightgray")      # sovereign debt crisis: Aug 2011 - Feb 2013
rect(date[134], u[3], date[135], u[4], border = 0, col = "lightgray")      # legislative change: Feb 2015
rect(date[194], u[3], date[195], u[4], border = 0, col = "lightgray")      # covid: March 2020 -
lines(date, rho_const_ts[-len_m], lwd=2, col ="blue")
lines(date, rho_CS_ts[-len_m], lwd=2, col ="red")
lines(date, rho_BF_ts[-1], lwd=2)
#legend("bottomright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))
box()


# plot gamma

#plot(date, gamma_CS_ts[-len_m], lwd=2, main=expression(hat(gamma)[t]), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red", ylim=c(-1,3))
plot(date, gamma_CS_ts[-len_m], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red", ylim=c(-1,3))
axis(side = 1, at = date[ix], labels = labs[ix], tcl = -0.7, cex.axis = 0.7)
u <- par("usr")
#plot green rect - notice that x-coordinates are defined by date points
rect(date[49], u[3], date[64], u[4], border = 0, col = "lightgray")      # global financial crisis
#rect(date[92], u[3], date[110], u[4], border = 0, col = "lightgray")      # sovereign debt crisis: Aug 2011 - Feb 2013
rect(date[134], u[3], date[135], u[4], border = 0, col = "lightgray")      # legislative change: Feb 2015
#rect(date[194], u[3], date[195], u[4], border = 0, col = "lightgray")      # covid: March 2020 -
lines(date, gamma_const_ts[-len_m], lwd=2, col ="blue")
lines(date, gamma_CS_ts[-len_m], lwd=2, col ="red")
lines(date, gamma_BF_ts[-1], lwd=2)
#legend("bottomright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))
box()


# plot observed TS

date_full <- seq(as.Date("2004/1/1"), by = "month", length.out = len_m)
fmt_full <- "%Y"
ix_full <- seq(1, length(date_full), 12)
labs_full <- format(date_full, fmt_full)

#plot(date_full, waves[,1], lwd=1, main="Observed time series", ylab=" ", xlab=" ", type="l", xaxt ="n", col ="black", ylim=c(1e+05, 8e+05))
plot(date_full, waves[,1], lwd=1, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="black", ylim=c(1e+05, 8e+05))
axis(side = 1, at = date[ix_full], labels = labs[ix_full], tcl = -0.7, cex.axis = 0.7)
lines(date_full, waves[,2], lwd=1)
lines(date_full, waves[,3], lwd=1)
lines(date_full, waves[,4], lwd=1)
lines(date_full, waves[,5], lwd=1)
lines(date_full, CC, lwd=1, col="red")
#legend("bottomright", legend = c("Labour force survey", "Claimant counts"),col = c("black","red"), lty = c(1, 1), lwd=c(1,1))



# plot estimated variances

date_cut <- date[-c(1:(30-13))]

ix_cut <- seq(1, length(date_cut), 12)
labs_cut <- format(date_cut, fmt)

# theta
#plot(date_cut, P_theta_biv_splines[-length(P_theta_biv_splines)], lwd=2, main=expression(hat(var)(hat(theta)[yt])), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
plot(date_cut, P_theta_biv_splines[-length(P_theta_biv_splines)], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
axis(side = 1, at = date_cut[ix], labels = labs_cut[ix], tcl = -0.7, cex.axis = 0.7)
lines(date_cut, P_theta_biv[-length(P_theta_biv)], lwd=2, col ="blue")
lines(date_cut, P_theta_biv_BF, lwd=2)
#legend("topright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))

# L
#plot(date_cut, P_L_biv_splines[-length(P_L_biv_splines)], lwd=2, main=expression(hat(var)(hat(L)[yt])), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
plot(date_cut, P_L_biv_splines[-length(P_L_biv_splines)], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
axis(side = 1, at = date_cut[ix], labels = labs_cut[ix], tcl = -0.7, cex.axis = 0.7)
lines(date_cut, P_L_biv[-length(P_L_biv)], lwd=2, col ="blue")
lines(date_cut, P_L_biv_BF, lwd=2)
#legend("bottomright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))

# R
#plot(date_cut, P_R_biv_splines[-length(P_R_biv_splines)], lwd=2, main=expression(hat(var)(hat(R)[yt])), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
plot(date_cut, P_R_biv_splines[-length(P_R_biv_splines)], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
axis(side = 1, at = date_cut[ix], labels = labs_cut[ix], tcl = -0.7, cex.axis = 0.7)
lines(date_cut, P_R_biv[-length(P_R_biv)], lwd=2, col ="blue")
lines(date_cut, P_R_biv_BF, lwd=2)
#legend("topright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))



# plot estimated state variables

# theta
#plot(date_cut, theta_biv_splines[-length(theta_biv_splines)], lwd=2, main=expression(hat(theta)[yt]), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
plot(date_cut, theta_biv_splines[-length(theta_biv_splines)], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
axis(side = 1, at = date_cut[ix], labels = labs_cut[ix], tcl = -0.7, cex.axis = 0.7)
lines(date_cut, theta_biv[-length(theta_biv)], lwd=2, col ="blue")
lines(date_cut, theta_biv_BF, lwd=2)
#legend("topright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))

# L
#plot(date_cut, L_biv_splines[-length(L_biv_splines)], lwd=2, main=expression(hat(L)[yt]), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
plot(date_cut, L_biv_splines[-length(L_biv_splines)], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
axis(side = 1, at = date_cut[ix], labels = labs_cut[ix], tcl = -0.7, cex.axis = 0.7)
lines(date_cut, L_biv[-length(L_biv)], lwd=2, col ="blue")
lines(date_cut, L_biv_BF, lwd=2)
#legend("bottomright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))

# R
#plot(date_cut, R_biv_splines[-length(R_biv_splines)], lwd=2, main=expression(hat(R)[yt]), ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
plot(date_cut, R_biv_splines[-length(R_biv_splines)], lwd=2, ylab=" ", xlab=" ", type="l", xaxt ="n", col ="red")
axis(side = 1, at = date_cut[ix], labels = labs_cut[ix], tcl = -0.7, cex.axis = 0.7)
lines(date_cut, R_biv[-length(R_biv)], lwd=2, col ="blue")
lines(date_cut, R_biv_BF, lwd=2)
#legend("topright", legend = c("Constant","Cubic splines", "Bootstrap filter"),col = c("blue","red","black"), lty = c(1, 1,1), lwd=c(2,2,2))



#### Comparison of methods ####

mean(P_L_biv_splines)/mean(P_L_biv)
mean(P_L_biv_BF)/mean(P_L_biv)

mean(P_R_biv_splines)/mean(P_R_biv)
mean(P_R_biv_BF)/mean(P_R_biv)

mean(P_theta_biv_splines)/mean(P_theta_biv)
mean(P_theta_biv_BF)/mean(P_theta_biv)



#### residuals diagnostics ####

eta <- res_BF[,(31-13):(len_m-1)]

# serial correlation
p_values_dt <- matrix(NA,2,6)

half <- trunc(ncol(eta)/2)

# univariate normality
for (i in 1:6){
  p_values_dt[1,i] <- as.numeric(shapiro.test(eta[i,])[2])      # normality
  F_cv <- sum((eta[i,(ncol(eta)-half+1):ncol(eta)])^2)/sum((eta[i,1:half])^2)
  p_values_dt[2,i] <- pf(F_cv, half, half, lower.tail = F, log.p = FALSE)      # heteroskedasticity
}




