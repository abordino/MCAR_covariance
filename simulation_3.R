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
n = 200
M = 50

little_power = c()
our_power = c()
R = c()

# Select the copula
cp = claytonCopula(param = c(1), dim = 5)

# Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
P = mvdc(copula = cp, margins = c("exp", "exp", "exp", "exp", "exp"),
         paramMargins = list(list(1), list(1), list(1), list(1), list(1)))

X = rMvdc(3*n, P)

for(p in seq(0.03, 0.3, length.out=10)){
  
  
  X = delete_MCAR(X, p, c(1,4,5))
  
  R = c(R, computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS)$R)
  
  
  ###### SAMPLE LEVEL, REPEATING THE TEST M TIMES #######
  little_decisions = c()
  our_decisions = c()
  for (i in 1:100){
    X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
    X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
    X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
    
    # create dataset with NA's
    columns = c("X1","X2","X3") 
    X = data.frame(matrix(nrow = 3*n, ncol = length(columns)))
    X[1:n, c("X1", "X2")] = X1
    X[(n+1):(2*n), c("X2", "X3")] = X2
    X[(2*n+1):(3*n), c("X1", "X3")] = X3
    
    little_decisions = c(little_decisions, little_test(X, alpha))
    our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
  }
  
  little_power = c(little_power, mean(little_decisions))
  our_power = c(our_power, mean(our_decisions))
}

plot(R, little_power, col="green", ylim = c(0,1), pch=18)
points(R, our_power, col="blue", pch=19)
abline(h = alpha, col="red")