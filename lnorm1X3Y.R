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
d = 8

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with lnormal and t marginals
P = mvdc(copula = cp, margins = c(rep("lnorm",d)),
         paramMargins = rep(list(1),d) )
data = rMvdc(n, P)

little_power = function(p){
  little_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MCAR(data, p, c(1,2,3,4,5))
    little_decision[i] = little_test(X, alpha = 0.05)
  }
  
  return(mean(little_decision, na.rm=TRUE))
}

bootstrap_power = function(p){
  our_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MCAR(data, p, c(1,2,3,4,5))
    our_decision[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  
  return(mean(our_decision, na.rm=TRUE))
}


######## USING FOREACH AND DORNG
registerDoFuture()
plan(multicore)
RNGkind("L'Ecuyer-CMRG")
set.seed(232)

start.time = Sys.time()

xxx = seq(0.03, 0.3, length.out=7)

little_power = foreach(t1 = xxx, .combine = 'c') %dorng% little_power(t1)
our_power = foreach(t1 = xxx, .combine = 'c') %dorng% bootstrap_power(t1)

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("Pictures/MCARlnorm_3X5Y.png")
plot(seq(0.03, 0.3, length.out=7), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's type I error", "Our type I error"),
       col = c("green", "blue"),
       pch = c(18, 18))
dev.off()


######### MAR ##############################
alpha = 0.05
n = 200
MC = 100
d = 8

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with lnormal and t marginals
P = mvdc(copula = cp, margins = c(rep("lnorm",d)),
         paramMargins = rep(list(1),d) )
data = rMvdc(n, P)

little_power = function(p){
  little_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MAR_1_to_x(data, p, c(1,2,3,4,5), cols_ctrl = c(6,6,7,8,8), x=9)
    little_decision[i] = little_test(X, alpha = 0.05)
  }
  
  return(mean(little_decision))
}

bootstrap_power = function(p){
  our_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MAR_1_to_x(data, p, c(1,2,3,4,5), cols_ctrl = c(6,6,7,8,8), x=9)
    our_decision[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  
  return(mean(our_decision))
}


######## USING FOREACH AND DORNG
registerDoFuture()
plan(multicore)
RNGkind("L'Ecuyer-CMRG")
set.seed(232)

start.time = Sys.time()

xxx = seq(0.03, 0.3, length.out=7)

little_power = foreach(t1 = xxx, .combine = 'c') %dorng% little_power(t1)
our_power = foreach(t1 = xxx, .combine = 'c') %dorng% bootstrap_power(t1)

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("Pictures/MAR_lnorm_3X5Y.png")
plot(seq(0.03, 0.3, length.out=7), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("bottomright",
       legend = c("Little's power", "Our power"),
       col = c("green", "blue"),
       pch = c(18, 18))
dev.off()

######### MAR ##############################
alpha = 0.05
n = 200
MC = 100
d = 8

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with lnormal and t marginals
P = mvdc(copula = cp, margins = c(rep("lnorm",d)),
         paramMargins = rep(list(1),d) )
data = rMvdc(n, P)

little_power = function(p){
  little_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MAR_rank(data, p, c(1,2,3,4,5), cols_ctrl = c(6,6,7,8,8))
    little_decision[i] = little_test(X, alpha = 0.05)
  }
  
  return(mean(little_decision))
}

bootstrap_power = function(p){
  our_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MAR_rank(data, p, c(1,2,3,4,5), cols_ctrl = c(6,6,7,8,8))
    our_decision[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  
  return(mean(our_decision))
}


######## USING FOREACH AND DORNG
registerDoFuture()
plan(multicore)
RNGkind("L'Ecuyer-CMRG")
set.seed(232)

start.time = Sys.time()

xxx = seq(0.03, 0.3, length.out=7)

little_power = foreach(t1 = xxx, .combine = 'c') %dorng% little_power(t1)
our_power = foreach(t1 = xxx, .combine = 'c') %dorng% bootstrap_power(t1)

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("Pictures/MAR_lnorm_rank_3X5Y.png")
plot(seq(0.03, 0.3, length.out=7), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("bottomright",
       legend = c("Little's power", "Our power"),
       col = c("green", "blue"),
       pch = c(18, 18))
dev.off()
