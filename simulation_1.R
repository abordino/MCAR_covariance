setwd("~/Documents/phd/MCAR_covariance/simulations")
source("simul_general.R")
source("little_test.R")
source("bootstrap_test.R")
library(naniar)
library(missMethods)
library(MASS)
library(norm)

######### example 1 on 3-cycle ############
alpha = 0.05
n = 200
M = 1000
t3 = pi/6
t2 = pi/4

little_power = c()
our_power = c()
R = c()

for(t1 in seq(t2-t3, pi - t3, length.out=20)){

  #### POPULATION LEVEL ######
  SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
  for(j in 1:3){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }
  
  SigmaS[[1]][1,2] = cos(t1)
  SigmaS[[1]][2,1] = cos(t1)
  SigmaS[[2]][1,2] = cos(t2)
  SigmaS[[2]][2,1] = cos(t2)
  SigmaS[[3]][1,2] = cos(t3)
  SigmaS[[3]][2,1] = cos(t3)
  
  R = c(R, computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS)$R)
  
  
  ###### SAMPLE LEVEL, REPEATING THE TEST M TIMES ########
  little_decisions = c()
  our_decisions = c()
  
  for(m in 1:M){
    print(m)
    X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
    X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
    X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
    
    # create dataset with NA's
    columns = c("X1","X2","X3") 
    X = data.frame(matrix(nrow = 3*n, ncol = length(columns)))
    X[1:n, c("X1", "X2")] = X1
    X[(n+1):(2*n), c("X2", "X3")] = X2
    X[(2*n+1):(3*n), c("X1", "X3")] = X3
    
    # run the two tests
    little_decisions = c(little_decisions, little_test(X, alpha))
    our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
  }
  
  little_power = c(little_power, mean(little_decisions))
  our_power = c(our_power, mean(our_decisions))
}

plot(R, little_power, col="green", ylim = c(0,1), pch=18)
points(R, our_power, col="blue", pch=19)
abline(h = alpha, col="red")

