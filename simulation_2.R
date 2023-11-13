setwd("~/Documents/phd/MCAR_covariance/simulations")
source("simul_general.R")
source("little_test.R")
source("bootstrap_test.R")
library(naniar)
library(missMethods)
library(MASS)
library(norm)
library(foreach)
library(doParallel)
library(copula)
library(missMethods)
library(misty)
library(dplyr)

######### example 1 on 3-cycle ############
alpha = 0.05
n = 1000
M = 1000

little_errorI = c()
our_errorI = c()

d = 6
A = matrix(runif(d^2)*2-1, ncol=d) 
Sigma = cov2cor(t(A) %*% A)
data = mvrnorm(n, rep(0,d), Sigma)

for(p in seq(0.03, 0.3, length.out=10)){
  
  little_decisions = c()
  our_decisions = c()
  
  for (i in 1:M){
    X = delete_MCAR(data, p, c(1,5,6))
    print("------------------------------------------------------------------")
    print(i)
    print("------------------------------------------------------------------")
    
    # if something is singular, skip iteration
    skip_to_next = FALSE
    try({
      little_decisions = c(little_decisions, little_test(X, alpha))
      our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
    }, 
    silent = TRUE)
  
  little_errorI = c(little_errorI, mean(little_decisions))
  our_errorI = c(our_errorI, mean(our_decisions))
  }
}

plot(seq(0.03, 0.3, length.out=10), our_errorI, col="blue", ylim = c(0,1), pch=18)
points(seq(0.03, 0.3, length.out=10), little_errorI, col="green", pch=19)
abline(h = alpha, col="red")