setwd("~/Documents/phd/MCAR_covariance/MCAR/")
source("computeR.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)

######### 3-cycle: setting 1 ############
alpha = 0.05
n = 200
M = 100

# definition of the block 3-cycle
d = 3

# definition of P
P = 0.5*matrix(runif((d)^2)*2-1, ncol=d)

while(norm(P, type="2") > 0.75){
  P = 0.5*matrix(runif((d)^2)*2-1, ncol=d)
}

normP = norm(P, type="2")

little_power = c()
little_power_cov = c()
our_power = c()
R = c()

for(beta in seq(1-2*normP^2, 0.97, length.out = 10)){
  
  #### POPULATION LEVEL ######
  SigmaS=list()
  
  SigmaS[[1]]=diag(2*d)
  SigmaS[[1]][1:d, (d+1):(2*d)] = P
  SigmaS[[1]][(d+1):(2*d), 1:d] = t(P)
  
  SigmaS[[2]]=diag(2*d)
  SigmaS[[2]][1:d, (d+1):(2*d)] = -P
  SigmaS[[2]][(d+1):(2*d), 1:d] = -t(P)
  
  SigmaS[[3]]=diag(2*d)
  SigmaS[[3]][1:d, (d+1):(2*d)] = beta*diag(d)
  SigmaS[[3]][(d+1):(2*d), 1:d] = beta*diag(d)
  
  result = computeR(list(c(1:(2*d)),c(1:d,(2*d+1):(3*d)),c((d+1):(3*d))),
                    SigmaS = SigmaS)
  
  R = c(R, result$R)
  
  
  ###### SAMPLE LEVEL, REPEATING THE TEST M TIMES #######
  little_decisions = c()
  little_decisions_cov = c()
  our_decisions = c()
  for (i in 1:M){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X1 = mvrnorm(n, rep(0, 2*d), SigmaS[[1]])
    X2 = mvrnorm(n, rep(0, 2*d), SigmaS[[2]])
    X3 = mvrnorm(n, rep(0, 2*d), SigmaS[[3]])
    
    X = data.frame(matrix(nrow = 3*n, ncol = length(columns)))
    X[1:n, c(1:(2*d))] = X1
    X[(n+1):(2*n), c(1:d,(2*d+1):(3*d))] = X2
    X[(2*n+1):(3*n), c((d+1):(3*d))] = X3
    X = as.matrix(X)
    
    ### run little's test
    little_decisions = c(little_decisions, little_test(X, alpha))
    little_decisions_cov = c(little_decisions_cov, little_test(X, alpha, "cov"))
    
    ### run our tests
    our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
  }
  
  little_power = c(little_power, mean(little_decisions))
  little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  our_power = c(our_power, mean(our_decisions))
}

plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "")
points(R, little_power_cov, col="orange", pch=18)
points(R, our_power, col="blue", pch=18)
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))