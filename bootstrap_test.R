setwd("~/Documents/phd/MCAR_covariance/simulations")
source("simul_general.R")
library(naniar)
library(missMethods)
library(MASS)

######### example 1 on 3-cycle ############
alpha = 0.05
n = 200
M = 1000
t2 = pi/6


little_power = c()
our_power = c()
for(t1 in seq(t2, pi/2, length.out=10)){
  
  little_decisions = c()
  our_decisions = c()
  
  SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
  for(j in 1:3){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }
  
  SigmaS[[1]][1,2] = cos(t1)
  SigmaS[[1]][2,1] = cos(t1)
  SigmaS[[2]][1,2] = cos(t2)
  SigmaS[[2]][2,1] = cos(t2)
  SigmaS[[3]][1,2] = 1
  SigmaS[[3]][2,1] = 1
  
  for(m in 1:M){
    print(m)
    X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
    X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
    X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
    
    #### compute Little's test statistics
    columns = c("X1","X2","X3") 
    X = data.frame(matrix(nrow = 3*n, ncol = length(columns)))
    X[1:n, c("X1", "X2")] = X1
    X[(n+1):(2*n), c("X2", "X3")] = X2
    X[(2*n+1):(3*n), c("X1", "X3")] = X3
    
    little_d = (mcar_test(X)$p.value < alpha)
    little_decisions = c(little_decisions, little_d)
    
    ##### compute our bootstrap test
    SigmaShat = list()
    SigmaShat[[1]] = cor(X1)
    SigmaShat[[2]] = cor(X2)
    SigmaShat[[3]] = cor(X3)
    
    # find the closest cov matrix to SigmaS
    result = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaShat)
    Q_hat = result$Sigma/(1-result$R)
    QS_hat = list()
    QS_hat[[1]] = Q_hat[c(1,2),c(1,2)]
    QS_hat[[2]] = Q_hat[c(2,3),c(2,3)]
    QS_hat[[3]] = Q_hat[c(1,3),c(1,3)]
    
    R_hat = computeR(list(c(1,2),c(2,3), c(1,3)), QS_hat)$R
    
    #### create bootstrap samples from X
    B = 99
    sum_indicator = 0
    for (b in 1:B){
      SigmaShat_b = list()
      SigmaShat_b[[1]] = cor(X1[sample(1:n, n, replace = T),])
      SigmaShat_b[[2]] = cor(X2[sample(1:n, n, replace = T),])
      SigmaShat_b[[3]] = cor(X3[sample(1:n, n, replace = T),])
      
      # find the closest cov matrix to SigmaS_b
      result = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaShat_b)
      Q_hat_b = result$Sigma/(1-result$R)
      QS_hat_b = list()
      QS_hat_b[[1]] = Q_hat_b[c(1,2),c(1,2)]
      QS_hat_b[[2]] = Q_hat_b[c(2,3),c(2,3)]
      QS_hat_b[[3]] = Q_hat_b[c(1,3),c(1,3)]
      
      R_hat_b = computeR(list(c(1,2),c(2,3), c(1,3)), QS_hat_b)$R
      sum_indicator = sum_indicator + as.numeric(R_hat_b <= R_hat)
    }
    
    our_d = (1 + sum_indicator <= alpha*(B+1))
    our_decisions = c(our_decisions, our_d)
  }
  little_power = c(little_power, mean(little_decisions))
  our_power = c(our_power, mean(our_decisions))
}

R = (cos(t2) - cos(seq(t2, pi/2, length.out=10)))/2
plot(R, our_power, col="orange", ylim = c(0,1), pch=18)
par(new=TRUE)
plot(R, little_power, col="green", ylim = c(0,1), pch=18)
abline(h = alpha, col="red")

