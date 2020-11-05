## Packages needed to run the script ##

library(ucminf)
library(magic)     
library(MASS)           
library(Rcpp)
library(RcppArmadillo)
library(zoo)
library(optimr)



## Functions ##


# Create the cubic splines weights
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


# Link function that bounds its argument between -1 and 1
link <- function(x){
  y <- x/sqrt(1+x^2)
  return(y)
}


# Kalman filter estimation with known values for the time-varying correlation
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
  
  logl <- 0
  
  for (i in 1:len){ 
    
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
    if ((NaN %in% Fmatrix)==T){
      logl<- -P10[1]
    } else {
      svdFmatrix <- svd(Fmatrix)
      Fmatrix_inv <- svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)
      if (opti==F){
        Bmatrix <- chol(Fmatrix_inv)
        st_for[,i] <- Bmatrix%*%epshatoutofsample
      }
      Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv 
      xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%epshatoutofsample 
      epshatinsample <- y[,i]-Z%*%xtt[,i] 
      Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%Z%*%Pttm1[[i]] 
      
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
      
      if (outofsample) {
        if (i <= d ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > d ){ 
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix_inv%*%epshatoutofsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      } else {
        if (i <= d ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > d ){ 
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


# Kalman filter estimation of the cubic splines model
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
  
  logl <- 0
  
  for (i in 1:len){ 
    
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
    if ((NaN %in% Fmatrix)==T){
      logl<- -P10[1]
    } else {
      svdFmatrix <- svd(Fmatrix)
      Fmatrix_inv <- svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)
      if (opti==F){
        Bmatrix <- chol(Fmatrix_inv)
        st_for[,i] <- Bmatrix%*%epshatoutofsample
      }
      Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv 
      xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%epshatoutofsample 
      epshatinsample <- y[,i]-Z%*%xtt[,i] 
      Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%Z%*%Pttm1[[i]] 
      
      if (hyper_tan==T){  
        if (restricted==T){
          R[32,2] <- tanh(par[length(par)-k+1])      
          R[2,32] <- tanh(par[length(par)-k+1]) 
        } else {
          R[32,2] <- tanh(W[i,1:k]%*%par[(length(par)-k+1):length(par)])      
          R[2,32] <- tanh(W[i,1:k]%*%par[(length(par)-k+1):length(par)]) 
        }
      } else {
        if (restricted==T){
          R[32,2] <- link(par[length(par)-k+1])      
          R[2,32] <- link(par[length(par)-k+1]) 
        } else {
          R[32,2] <- link(W[i,1:k]%*%par[(length(par)-k+1):length(par)])      
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
        } else if (i > d ){ 
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix_inv%*%epshatoutofsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      } else {
        if (i <= d ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > d ){ 
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


# Kalman filter estimation of the model with time-constant correlation
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
  
  logl <- 0
  
  for (i in 1:len){ 

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
    if ((NaN %in% Fmatrix)==T){
      logl<- -P10[1]
    } else {
      svdFmatrix <- svd(Fmatrix)
      Fmatrix_inv <- svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)
      if (opti==F){
        Bmatrix <- chol(Fmatrix_inv)
        st_for[,i] <- Bmatrix%*%epshatoutofsample
      }
      Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv 
      xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%epshatoutofsample 
      epshatinsample <- y[,i]-Z%*%xtt[,i] 
      Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix_inv%*%Z%*%Pttm1[[i]] 
      Pttm1[[i+1]] <- Tmatrix%*%Pttm1[[i]]%*%t(Tmatrix-Kg%*%Z)+Q 
      xttm1[,i+1] <- Tmatrix%*%xttm1[,i] + Kg%*%epshatoutofsample
      
      if (outofsample) {
        if (i <= (30-13) ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > (30-13) ){ 
          logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix_inv%*%epshatoutofsample
          if ((NaN %in% logl)==T){
            logl<- -P10[1]
          }
        }
      } else {
        if (i <= (30-13) ){
          logl <- logl - nrow(y)/2*log(2*pi)
        } else if (i > (30-13) ){ 
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


# Estimation of the static parameters of the model by indirect inference
IndInf_CC <- function(par, sim_h, beta_hat, len, eps_sim, VARgamma, opti,outofsample,parP10,nstates,d,k,W,restricted,hyper_tan,se,states_noerr,nvar,init_val_CC){ 
  sigma_Ry <- par[1]
  sigma_omegay <- par[2]
  sigma_lambda <- par[3]
  sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
  sigma_Rx <- par[9]
  sigma_omegax <- par[10]  
  betatilde <- matrix(NA, nrow=sim_h, ncol=length(par))
  
  possibleError <- tryCatch(for (h in 1:sim_h){
    
    # generate time-varying corr according to a random walk model
    gammaII <- cumsum(exp(par[length(par)])*eps_sim[2,(h*len-len+1):(h*len)])
    if (hyper_tan == T){
      rhoII <- tanh(gammaII)
    } else {
      rhoII <- link(gammaII)
    }
    
    # generate state innovations
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
    
    # generate innovations in the observation equation
    nu <- exp(par[11])*eps_sim[1,(h*len-len+1):(h*len)]  
    Msel <- diag(1,nvar)      # selection matrix of the innovations in the observation equation
    Msel <- Msel[,-c(1:5)]
    
    # generate the observable variables
    y <- matrix(NA,nvar,len)
    for (i in 1:len){
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

  return(opt.fun)
}


                
## Upload cpp script ##

Cpp_IndInf_ssm_levels <- sourceCpp("estimation_CC_IndInf_BF.cpp")

                
                
## Estimation of models ##

# Initial values of static parameter
init_val_CC <- c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                 log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),
                 log(3000),log(0.02),0,log(1000), 0.21)



# Time-constant correlation #

objopt_CC_const <- ucminf(par=init_val_CC, KF_CC_const, y=y, se=se, opti=T, outofsample=T, parP10=1000000000000, nstates=43, 
                          hyper_tan=hyper_tan, hessian=2, control=list(gradstep = c(1e-2, 1e-3)))
                
par_biv_const <- objopt_CC_const$par

obj_const <- KF_CC_const(par_biv_const,y=y,se=se,opti=F,outofsample=T,parP10=1000000000000,nstates=43, hyper_tan=hyper_tan)
                
# estimated time-constant correlation
gamma_hat_const <- par_biv_const[11]
if (hyper_tan==T){
  rho_hat_const <- tanh(gamma_hat_const)
} else {
  rho_hat_const <- link(gamma_hat_const)
}


# Cubic splines model #

AIC <- rep(NA,2)
BIC <- rep(NA,2)

for (choice_k in 1:2){
  
  if (choice_k==1){ 
    knots <- as.numeric(quantile(c(1:len_m), probs = seq(0, 1, 0.25)))
  } else {
    knots <- as.numeric(quantile(c(1:len_m), probs = seq(0, 1, 1/7)))
  }
  
  W <- Weights(len=len_m,knots=knots)$W
  
  objopt_CC_splines <- ucminf_rcpp_splines(init_val=c(init_val_CC[-11], rep(0,length(knots))), y=y, se=as.matrix(se), opti=T, outofsample=T,
                                           parP10=1000000000000, nstates=43, d=30-13, k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, 
                                           control=list(gradstep = c(1e-2, 1e-3)))  
  
  AIC[choice_k] <- 1/len_m*(2*objopt_CC_splines$value + 2*(30-13 + length(c(init_val_CC[-11],rep(0,length(knots))))))
  BIC[choice_k] <- 1/len_m*(2*objopt_CC_splines$value + log(len_m)*(30-13 + length(c(init_val_CC[-11],rep(0,length(knots))))))
  
}

if (which.min(BIC) == 1){
  knots <- as.numeric(quantile(c(1:len_m), probs = seq(0, 1, 0.25)))
  
  W <- Weights(len=len_m,knots=knots)$W
  
  objopt_CC_splines <- ucminf_rcpp_splines(init_val=c(init_val_CC[-11],rep(0,length(knots))), y=y, se=as.matrix(se), opti=T, outofsample=T,
                                           parP10=1000000000000, nstates=43, d=30-13, k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, 
                                           control=list(gradstep = c(1e-2, 1e-3)))
}
                
par_biv_splines <- objopt_CC_splines$par

obj_splines <- KF_CC_splines(objopt_CC_splines$par, y=y, se=se, opti=F, outofsample=T, parP10=1000000000000, nstates=43,
                             k=length(knots), restricted=F, W=W, d=30-13, hyper_tan=hyper_tan)

# estimated time-varying correlation
gamma_hat_splines <- W[,1:length(knots)]%*%par_biv_splines[(length(par_biv_splines)-length(knots)+1):length(par_biv_splines)]
if (hyper_tan==T){
  rho_hat_splines <- tanh(gamma_hat_splines)
} else {
  rho_hat_splines <- link(gamma_hat_splines)
}


# Indirect inference #

# state variables that do not have an error term in the transition equation
states_noerr <- c(1,23:31)
Rsel <- diag(1,43)      # selection matrix of the innovations in the transition equation
Rsel <- Rsel[,-states_noerr]
Msel <- diag(1,6)      # selection matrix of the innovations in the observation equation
Msel <- Msel[,-c(1:5)]

sim_h <- 5
set.seed(2609)
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


IndInf_result <- optimr(par=c(par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))], log(grid_sigma_gamma[which.min(obj_grid)])), IndInf_CC, 
                        sim_h=sim_h, beta_hat=par_biv_splines[-c((length(par_biv_splines)-length(knots)+1):length(par_biv_splines))], len=len_m, eps_sim=eps_sim, VARgamma = var(par_biv_splines[13:length(par_biv_splines)]),
                        opti=T, outofsample=T, parP10=1000000000000, nstates=43, d=30-13, k=length(knots), W=W, hyper_tan=hyper_tan, restricted=F, se=se, states_noerr=states_noerr,nvar=6,init_val_CC=init_val_CC,
                        control=list(trace=T), method="CG")


# Bootstrap filter #

draw_m <- 5000
set.seed(1006)
                
results_BF <- boot_filter_CC_rcpp(draw_m, IndInf_result$par, y, as.matrix(se), nstates=43, hyper_tan, Rsel, states_noerr, init_gamma = par_biv_const[11])

KF_final <- KF_CC_known_corr(par=IndInf_result$par[1:12],y=y[,-len_m],se,opti=F,outofsample=T,parP10=1000000000000,nstates=43,d=30-13,
                             gamma_draw = results_BF$att_BF[2:len_m], hyper_tan = hyper_tan)
                
# estimated time-varying correlation
if (hyper_tan == T){
  rho_hat_BF = tanh(results_BF$att_BF)
} else {
  rho_hat_BF = link(results_BF$att_BF)
}
