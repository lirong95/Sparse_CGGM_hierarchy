
## soft-threshold function
S_soft <- function(z, lambda){
  return((abs(z) - lambda)*(abs(z) - lambda > 0)*sign(z))
}

## make the matrix to be symmetric
Symmetrize <- function(X){
  p = dim(X)[1]
  for(i in 1:p){
    for(j in i:p){
      if(X[i,j] < X[j, i]){
        X[j, i] = X[i, j]
      }else{
        X[i, j] = X[j, i]
      }
    }
  }
  return(X)
}

## negative log-likelihood function of L_2
L_Omega <- function(Sigma_yx, Sigma_xx, Sigma_yy, Omega_yy, Omega_yx){
  trace1 <- sum(diag(Sigma_yy %*% Omega_yy))
  trace2 <- sum(diag(t(Sigma_yx) %*% Omega_yx))
  trace3 <- sum(diag(Sigma_xx %*% t(Omega_yx) %*% ginv(Omega_yy) %*% Omega_yx))
  if(det(Omega_yy) > 0){
    return( -log(det(Omega_yy)) + trace1 + 2*trace2 + trace3)
  }else{
    return(Inf)
  } 
}

## update (Omega_YY, Omega_YX, S) in the iteration
#@lambda1: tuning parameter for the group effects
#@lambda2: tuning parameter for S
#@lambda3: tuning parameter for gamma (balance of sparsity and low-rank)
#@lambda4: tuning parameter for Omega_YX
#@maxiter: the maximum iteration of gradient descent algorithm
#@tol: the minimum difference between previous estimate and current estimate in gradient descent algorithm
#@zeta0: initial zeta in gradient descent algorithm
#@tau: Lagrangian parameter, default = 1
Update_Omega <- function(lambda1, lambda2, lambda3, lambda4, Sigma_yx, Sigma_xx, Sigma_yy,
                         R_prev, L_prev, Phi_prev,
                         Omega_yy_prev, Omega_yx_prev, S_prev,
                         maxiter, tol, zeta0, tau){
  p <- dim(Sigma_yx)[1]
  q <- dim(Sigma_yx)[2]
  tilde_Omega_yx <- Omega_yx_prev
  tilde_Omega_yy <- Omega_yy_prev
  tilde_S <- S_prev
  diff_value = 10
  t = 0
  zeta = zeta0
  diff_vec <- NULL
  while( diff_value >= tol  && t < maxiter ){
    tilde_Omega_yx_old <- tilde_Omega_yx 
    tilde_Omega_yy_old <- tilde_Omega_yy
    tilde_S_old <- tilde_S
    
    A.matrix <- 2*Sigma_yx + 2*ginv(tilde_Omega_yy_old) %*% tilde_Omega_yx_old %*% Sigma_xx
    B.matrix <- Sigma_yy - ginv(tilde_Omega_yy_old) - ginv(tilde_Omega_yy_old)%*%tilde_Omega_yx_old%*%Sigma_xx%*%t(tilde_Omega_yx_old)%*%ginv(tilde_Omega_yy_old)
    W_omega <- tilde_Omega_yy_old - zeta*B.matrix
    W_s <- tilde_S_old - zeta/tau*(tilde_S_old - R_prev - L_prev + tau*Phi_prev)
    tilde_Omega_yx_updated <- matrix(0, p, q)
    for(j in 1:p){
      for(m in 1:q){
        tilde_Omega_yx_updated[j,m] <- S_soft( (tilde_Omega_yx_old[j,m]-zeta*A.matrix[j,m]), zeta*lambda4)
      }
    }
    tilde_Omega_yy_updated <- tilde_S_updated <- matrix(0, p, p)
    for(j in 1:p){
      for(m in 1:p){
        u_bar <- norm(c( W_omega[j,m], S_soft(W_s[j,m], zeta*lambda2) ), type="2")
        group_effect <- max(0, (1 - zeta*lambda1/u_bar))
        tilde_Omega_yy_updated[j,m] <- group_effect * W_omega[j,m]
        tilde_S_updated[j,m] <- group_effect * S_soft(W_s[j,m], zeta*lambda2)
      }
    }
    L.like <- L_Omega(Sigma_yx, Sigma_xx, Sigma_yy, tilde_Omega_yy_updated, tilde_Omega_yx_updated)
    term.1 <- sum(diag(t(A.matrix) %*% (tilde_Omega_yx_updated - tilde_Omega_yx_old)))
    term.2 <- sum(diag(t(B.matrix) %*% (tilde_Omega_yy_updated - tilde_Omega_yy_old)))
    term.3 <- 1/(2*zeta)*(norm(tilde_Omega_yx_updated - tilde_Omega_yx_old, type="F")^2+norm(tilde_Omega_yy_updated - tilde_Omega_yy_old, type="F")^2)
    tilde.L.like <- L_Omega(Sigma_yx, Sigma_xx, Sigma_yy, tilde_Omega_yy_old, tilde_Omega_yx_old) + term.1 + term.2 + term.3
    if(L.like <= tilde.L.like){
      tilde_Omega_yx <- tilde_Omega_yx_updated
      for(j in 1:p){
        for(m in 1:p){
          if(j != m){
            tilde_Omega_yy[j,m] <- tilde_Omega_yy_updated[j,m]
            tilde_S[j,m] <- tilde_S_updated[j,m]
          }else{
            tilde_Omega_yy[j,m] <- tilde_Omega_yy_updated[j,m]
            tilde_S[j,m] <- tilde_S_updated[j,m]
          }
        }
      }
      tilde_Omega_yy <- Symmetrize(tilde_Omega_yy)
      tilde_S <- Symmetrize(tilde_S)
      tilde_Omega_yy[abs(tilde_Omega_yy) < 1e-4] = 0
      tilde_Omega_yx[abs(tilde_Omega_yx) < 1e-4] = 0
      tilde_S[abs(tilde_S) < 1e-4] = 0
      t = t + 1
      diff_value1 = norm(tilde_Omega_yx_old-tilde_Omega_yx,type="F")
      diff_value2 = norm(tilde_Omega_yy_old-tilde_Omega_yy,type="F")
      diff_value3 = norm(tilde_S_old-tilde_S,type="F")
      diff_value <- max(diff_value1, diff_value2, diff_value3)
      diff_vec <- c(diff_vec, diff_value) 
    }else{
      zeta <- 0.5*zeta
    }
  }
  res.list <- list()
  res.list$t = t; res.list$diff = diff_vec; res.list$zeta = zeta
  res.list$Omega_yx = tilde_Omega_yx
  res.list$Omega_yy = tilde_Omega_yy
  res.list$S = tilde_S
  return(res.list)
}

# update L in the iteration
Update_L <- function(lambda1, lambda2, lambda3, S_prev, R_prev, Phi_prev, tau){
  W_l <- S_prev - R_prev + tau * Phi_prev
  SaS <- W_l
  edecomp = eigen(SaS)
  Da = edecomp$values
  if(is.complex(Da)){SaS=Symmetrize(W_l); edecomp = eigen(SaS); Da = edecomp$values}
  Va = edecomp$vectors
  Da2 <- NULL
  for(j in 1:length(Da)){
    Da2[j] = max(0, Da[j] - tau*lambda3)
  }
  L = Va %*% diag(Da2) %*% t(Va)
  return(L)
}

# update R in the iteration
Update_R <- function(lambda1, lambda2, lambda3, Sigma_yy,
                     S_prev, L_prev, Phi_prev, tau){
  W_r <- S_prev - L_prev + tau * Phi_prev 
  Sa = tau * Sigma_yy - W_r
  edecomp = eigen(Sa)
  D = edecomp$values
  if(is.complex(D)){Sa=Symmetrize(Sa); edecomp = eigen(Sa); D = edecomp$values}
  V = edecomp$vectors
  D2 = (-D + sqrt(D^2 + 4*tau))/2
  R = V %*% diag(D2) %*% t(V)
  return(R)
}

## the iterative process to estimate
#@data: Y, the set of variables of interest
#@covarite: X, the set of additional variables
#@eps: the minimum difference between previous estimate and current estimate in ADMM
#@ niter: the maximum iteration of ADMM
#@initialize: the list of initial value
CLGGM <- function(data, covariate, lambda1, lambda2, lambda3, lambda4,
                  eps, niter, zeta0, tau, tol, maxiter,
                  initialization=TRUE, initialize){
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  q <- as.integer(dim(covariate)[2])
  Sigma_yy <- t(data) %*% data / n
  Sigma_xx <- t(covariate) %*% covariate / n
  Sigma_yx <- t(data) %*% covariate / n
  Omega_yy = initialize$Omega_yy
  Omega_yx = initialize$Omega_yx
  S = initialize$S
  L = matrix(0, p, p)
  R = S + L
  Phi = matrix(0, p, p)
  
  # start initialize
  t = 0
  diff_value = 10
  while( diff_value >= eps  && t < maxiter ){
    Omega_yy_old = Omega_yy
    Omega_yx_old = Omega_yx
    S_old = S
    L_old = L
    R_old = R
    Phi_old = Phi
    
    R <- Update_R(lambda1, lambda2, lambda3, Sigma_yy, S_old, L_old, Phi_old, tau)
    L <- Update_L(lambda1, lambda2, lambda3, S_old, R, Phi_old, tau)
    update.spagroup <- Update_Omega(lambda1, lambda2, lambda3, lambda4, Sigma_yx, Sigma_xx, Sigma_yy,
                                    R, L, Phi_old,
                                    Omega_yy_old, Omega_yx_old, S_old,
                                    maxiter, tol, zeta0, tau)
    Omega_yy <- update.spagroup$Omega_yy
    S <- update.spagroup$S
    Omega_yx <- update.spagroup$Omega_yx
    
    Phi <- Phi_old - (R - S + L)/tau
    
    t = t + 1
    diff_Omega_yy = norm(Omega_yy_old-Omega_yy,type="F")/(norm(Omega_yy_old, type="F")+0.0001)
    diff_S = norm(S_old-S,type="F")/(norm(S_old, type="F")+0.0001)
    diff_L = norm(L_old-L,type="F")/(norm(L_old, type="F")+0.0001)
    diff_R = norm(R_old-R,type="F")/(norm(R_old, type="F")+0.0001)
    diff_value <- max(diff_Omega_yy, diff_S, diff_L, diff_R)
  }
  Omega_yy <- Symmetrize(Omega_yy)
  S <- Symmetrize(S)
  
  Omega_yy[abs(Omega_yy) < 1e-2] = 0
  Omega_yx[abs(Omega_yx) < 1e-2] = 0
  S[abs(S) < 1e-2] = 0
  CLGGM_res <- list()
  CLGGM_res$Omega_yy <- Omega_yy; CLGGM_res$Omega_yx <- Omega_yx
  CLGGM_res$S <- S; CLGGM_res$L <- L
  CLGGM_res$t <- t
  return(CLGGM_res)
}

## BIC function
BIC <- function(S_est, L_est, Omega_yy_est, Omega_yx_est, data, covariate){
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  q <- as.integer(dim(covariate)[2])
  Sigma_yy <- t(data) %*% data / n
  Sigma_xx <- t(covariate) %*% covariate / n
  Sigma_yx <- t(data) %*% covariate / n
  nonz = sum(S_est!=0)+sum(Omega_yy_est!=0)+sum(Omega_yx_est!=0)
  if(det(S_est) > 0){
    l2.value <- -log(det(S_est-L_est)) + sum(diag(Sigma_yy%*%(S_est-L_est)))
  }else{ l2.value <- Inf}
  if(det(Omega_yy_est) > 0){
    l1.value <- -log(det(Omega_yy_est)) + sum(diag(Sigma_yy%*%Omega_yy_est)) + 2*sum(diag(t(Sigma_yx)%*%Omega_yx_est)) + sum(diag(Sigma_xx%*%t(Omega_yx_est)%*%solve(Omega_yy_est)%*%Omega_yx_est))
  }else{ l1.value <- Inf}
  bic <- 2*n*(l1.value + l2.value) + log(n)*nonz
  return(list(bic=bic, l1=l1.value, l2=l2.value, nonz=nonz, nonz1=sum(S_est!=0), nonz2=sum(Omega_yx_est!=0)))
}








