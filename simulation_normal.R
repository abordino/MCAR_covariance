setwd("~/Documents/phd/MCAR_covariance/simulations")
source("simul_general.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)
library(copula)
library(missMethods)

######### example 1 on 3-cycle ############
alpha = 0.05
n = 1000
M = 100

little_errorI = c()
our_errorI = c()

d = 5
A = matrix(runif(d^2)*2-1, ncol=d) 
Sigma = cov2cor(t(A) %*% A)
data = mvrnorm(n, rep(0,d), Sigma)


######## simulation of type I error, if data are MCAR 
for(p in seq(0.03, 0.3, length.out=10)){
  
  little_decisions = c()
  our_decisions = c()
  
  for (i in 1:M){
    X = delete_MCAR(data, p, c(1,2,5))
    print("------------------------------------------------------------------")
    print(i)
    print("------------------------------------------------------------------")
    
    little_decisions = c(little_decisions, little_test(X, alpha))
    our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
    
  }
  
  little_errorI = c(little_errorI, mean(little_decisions))
  our_errorI = c(our_errorI, mean(our_decisions))
}

plot(seq(0.03, 0.3, length.out=10), little_errorI, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "")
points(seq(0.03, 0.3, length.out=10), our_errorI, col="blue", pch=18)
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's type I error", "Our type I error"),
       col = c("green", "blue"),
       pch = c(18, 18))



######## simulation of power, if data are MAR with censoring 
little_power = c()
our_power = c()

d = 5
A = matrix(runif(d^2)*2-1, ncol=d) 
Sigma = cov2cor(t(A) %*% A)
data = mvrnorm(n, rep(0,d), Sigma)

for(p in seq(0.03, 0.3, length.out=10)){
  
  little_decisions = c()
  our_decisions = c()
  
  for (i in 1:M){
    X = delete_MAR_1_to_x(data, p, c(1,2,5), c(3,3,4), x = 9)
    print("------------------------------------------------------------------")
    print(i)
    print("------------------------------------------------------------------")
    
    little_decisions = c(little_decisions, little_test(X, alpha))
    our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
    
  }
  
  little_power = c(little_power, mean(little_decisions))
  our_power = c(our_power, mean(our_decisions))
}

plot(seq(0.03, 0.3, length.out=10), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "")
points(seq(0.03, 0.3, length.out=10), our_power, col="blue", pch=18)
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's power", "Our power"),
       col = c("green", "blue"),
       pch = c(18, 18))


