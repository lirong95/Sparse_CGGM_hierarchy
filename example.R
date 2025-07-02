library(mvtnorm)
library(MASS)

source("Sparse_CGGM_Hierarchy.R")

## Data generating process
set.seed(123456)
n <- 800
p <- 50
q <- 50
h <- 2
sub.numb <- 10
K_yy <- matrix(0, nrow=p, ncol=p)
K_yx <- matrix(0, nrow=p, ncol=q)
K_yeta <- matrix(0, nrow=p, ncol=h)
K_xx <- matrix(0, nrow=q, ncol=q)
K_etaeta <- matrix(0, nrow=h, ncol=h)
pp <- c(0, cumsum(rep(p/sub.numb, sub.numb)))
qq <- c(0, cumsum(rep(q/sub.numb, sub.numb)))
for(sub in 1:sub.numb){
  for(i in (pp[sub]+1):pp[sub+1]){
    for(j in (pp[sub]+1):pp[sub+1]){
      if(abs(i-j)==1){
        K_yy[i,j] = 0.4
      }
      if(i==j) K_yy[i,j] = 1
    }
  }
  for(i in (qq[sub]+1):qq[sub+1]){
    for(j in (qq[sub]+1):qq[sub+1]){
      if(abs(i-j)==1){
        K_xx[i,j] = 0.2
      }
      if(i==j) K_xx[i,j] = 1
    }
  }
  for(i in 1:h){
    for(j in 1:h){
      if(abs(i-j)==1){
        K_etaeta[i,j] = 0.2
      }
      if(i==j) K_etaeta[i,j] = 1
    }
  }
  all.comb <- NULL
  for(i in (pp[sub]+1):pp[sub+1]){
    for(j in (qq[sub]+1):qq[sub+1]){
      all.comb <- cbind(all.comb, c(i,j))
    }
  }
  ind.sample <- sample(1:ncol(all.comb), size=p/sub.numb)
  for(i in 1:length(ind.sample)){
    K_yx[all.comb[,ind.sample[i]][1], all.comb[,ind.sample[i]][2]] <- runif(min=c(-0.4,0.2), max=c(-0.2,0.4),1)
  }
  all.comb <- NULL
  for(i in (pp[sub]+1):pp[sub+1]){
    for(j in 1:h){
      all.comb <- cbind(all.comb, c(i,j))
    }
  }
  ind.sample <- sample(1:ncol(all.comb), size=5)
  for(i in 1:length(ind.sample)){
    K_yeta[all.comb[,ind.sample[i]][1], all.comb[,ind.sample[i]][2]] <- sample(c(0.1, -0.1), 1)
  }
  
}
K_xeta <- t(K_yx) %*% solve(K_yy) %*% K_yeta

S_star <- K_yy
marg.xeta <- rbind(cbind(K_xx, K_xeta), cbind(t(K_xeta), K_etaeta))
L_star <- cbind(K_yx, K_yeta) %*% solve(marg.xeta) %*% rbind(t(K_yx), t(K_yeta))
Omega_yy_star <- K_yy - K_yeta %*% solve(K_etaeta) %*% t(K_yeta)
Omega_yx_star <- K_yx - K_yeta %*% solve(K_etaeta) %*% t(K_xeta)
qr(L_star)$rank
Omega_yy_star[abs(Omega_yy_star) < 0.02] = 0
Omega_yx_star[abs(Omega_yx_star) < 0.02] = 0

K_whole <- rbind(cbind(K_yy, K_yx, K_yeta), 
                 cbind(t(K_yx), K_xx, K_xeta), 
                 cbind(t(K_yeta), t(K_xeta), K_etaeta))
set.seed(1234)
whole.data <- mvtnorm::rmvnorm(n, mean = rep(0, p+q+h), sigma = solve(K_whole))
Y <- whole.data[,1:p]
X <- whole.data[,(p+1):(p+q)]
eta <- whole.data[,(p+q+1):(p+q+h)]
data = Y
covariate = X
n <- dim(data)[1]
p <- dim(data)[2]
q <- dim(covariate)[2]
data.whole <- cbind(data, covariate)

## generate initial value
covmatrix <- cov(data.whole)
prematrix <- solve(covmatrix)
covdata <- cov(data)
initialize <- list()
initialize$Omega_yy <- prematrix[1:p, 1:p]
initialize$Omega_yy[abs(initialize$Omega_yy) < 0.01] = 0
initialize$Omega_yx <- prematrix[1:p, (p+1):(p+q)]
initialize$Omega_yx[abs(initialize$Omega_yx) < 0.01] = 0
initialize$S <- prematrix[1:p, 1:p]
initialize$S[abs(initialize$S) < 0.01] = 0

## estimation
res <- CLGGM(data, covariate, lambda1=0.05, lambda2=0.01, lambda3=0.02, lambda4=0.1,
             eps=5e-3, niter=40, zeta0=0.1, tau=1, tol=1e-2, maxiter=40,
             initialization=TRUE, initialize)




