setwd("~/Documents/GitHub/MCAR")
source("computeR.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
source("indexConsistency.R")
library(missMethods)
library(MASS)
library(norm)
library(foreach)
library(doSNOW)
library(doParallel)
library(future)
library(parallel)
library(foreach)
library(doRNG) 
library(doFuture)
library(future.apply)

######### MCAR ############
alpha = 0.05
n = 200
MC = 100
d = 3

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
P = mvdc(copula = cp, margins = c(rep("lnorm",d)),  #chisq #exp
         paramMargins = rep(list(1),d) )
data = rMvdc(n, P)


bootstrap_power = function(p){
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  our_decisions = logical(length = MC)
  for (i in 1:MC){

    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X = delete_MCAR(data, p, c(1,2))

    ### run our tests
    our_decisions[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  return(mean(our_decisions))
}

bootstrap_power_P = function(p){
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  our_decisions = logical(length = MC)
  for (i in 1:MC){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X = delete_MCAR(data, p, c(1,2))
    
    ### run our tests
    our_decisions[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "p")
  }
  return(mean(our_decisions))
}

little_power = function(p){
  little_decisions = logical(length = MC)
  for (i in 1:MC){

    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X = delete_MCAR(data, p, c(1,2))
    ### run little's tests
    little_decisions[i] = little_test(X, alpha = 0.05)
  }

  return(mean(little_decisions))
}

######## USING FOREACH AND DORNG
registerDoFuture()
plan(multicore)
RNGkind("L'Ecuyer-CMRG")
set.seed(232)

start.time = Sys.time()

xxx = seq(0.03, 0.3, length.out=7)

little_errorI = foreach(p = xxx, .combine = 'c') %dorng% little_power(p)
our_errorI = foreach(p = xxx, .combine = 'c') %dorng% bootstrap_power(p)
our_errorI_P = foreach(p = xxx, .combine = 'c') %dorng% bootstrap_power_P(p)

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("pictures/MCARlnorm_1X2Y.png")
plot(seq(0.03, 0.3, length.out=7), little_errorI, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_errorI, col="blue", pch=18, type = "b")
lines(seq(0.03, 0.3, length.out=7), our_errorI_P, col="cyan", pch=18, type = "b")
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's type I error", "Our type I error", "Our type I error P"),
       col = c("green", "blue", "cyan"),
       pch = c(18, 18, 18))
dev.off()

######### MAR ##############################
alpha = 0.05
n = 200
MC = 100
d = 3

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
P = mvdc(copula = cp, margins = c(rep("lnorm",d)),
         paramMargins = rep(list(1),d) )
data = rMvdc(n, P)


bootstrap_power = function(p){
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  our_decisions = logical(length = MC)
  for (i in 1:MC){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X = delete_MAR_1_to_x(data, p, c(1,2), cols_ctrl = c(3,3), x = 9)
    ### run our tests
    our_decisions[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  return(mean(our_decisions))
}

bootstrap_power_P = function(p){
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  our_decisions = logical(length = MC)
  for (i in 1:MC){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X = delete_MAR_1_to_x(data, p, c(1,2), cols_ctrl = c(3,3), x = 9)
    ### run our tests
    our_decisions[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "p")
  }
  return(mean(our_decisions))
}

little_power = function(p){
  little_decisions = logical(length = MC)
  for (i in 1:MC){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X = delete_MAR_1_to_x(data, p, c(1,2), cols_ctrl = c(3,3), x = 9)
    ### run little's tests
    little_decisions[i] = little_test(X, alpha = 0.05)
  }
  
  return(mean(little_decisions))
}

######## USING FOREACH AND DORNG
registerDoFuture()
plan(multicore)
RNGkind("L'Ecuyer-CMRG")
set.seed(232)

start.time = Sys.time()

xxx = seq(0.03, 0.3, length.out=7)

little_power = foreach(p = xxx, .combine = 'c') %dorng% little_power(p)
our_power = foreach(p = xxx, .combine = 'c') %dorng% bootstrap_power(p)
our_power_P = foreach(p = xxx, .combine = 'c') %dorng% bootstrap_power_P(p)


end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("pictures/MAR_lnorm_1X2Y.png")
plot(seq(0.03, 0.3, length.out=7), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_power, col="blue", pch=18, type = "b")
lines(seq(0.03, 0.3, length.out=7), our_power_P, col="cyan", pch=18, type = "b")
abline(h = alpha, col="red")
legend("bottomright",
       legend = c("Little's power", "Our power", "Our power P"),
       col = c("green", "blue", "cyan"),
       pch = c(18, 18, 18))
dev.off()