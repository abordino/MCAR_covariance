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
library(naniar)

######### MCAR ############
alpha = 0.05
n = 200
MC = 100
d = 5

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with chisqal and t marginals
P = mvdc(copula = cp, margins = c(rep("norm",d)),  #chisq #exp # norm
         paramMargins = rep(list(1),d) )

data = rMvdc(n, P)

xxx = seq(0.03, 0.3, length.out=7)

little_power = numeric(length = 7)
our_power = numeric(length = 7)
ind = 1
for (p in xxx){
  little_decision = logical(length = MC)
  our_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MCAR(data, p, c(1,2,5))
    
    little_decision[i] = mcar_test(data.frame(X))$p.value < alpha
    our_decision[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  
  little_power[ind] = mean(little_decision, na.rm = TRUE)
  our_power[ind] = mean(our_decision, na.rm = TRUE)
  ind = ind+1
}

plot(seq(0.03, 0.3, length.out=7), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's type I error", "Our type I error"),
       col = c("green", "blue"),
       pch = c(18, 18))


######### MAR ##############################
alpha = 0.05
n = 1000
MC = 100
d = 5

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with chisqal and t marginals
P = mvdc(copula = cp, margins = c(rep("norm",d)),
         paramMargins = rep(list(1),d) )

data = rMvdc(n, P)

xxx = seq(0.03, 0.3, length.out=7)

little_power = numeric(length = 7)
our_power = numeric(length = 7)
our_power_P = numeric(length = 7)
ind = 1
for (p in xxx){
  little_decision = logical(length = MC)
  our_decision = logical(length = MC)
  for (i in 1:MC){
    X = delete_MAR_1_to_x(data, p, c(1,2,5), cols_ctrl = c(3,4), x = 9)
    
    little_decision[i] = mcar_test(data.frame(X))$p.value < alpha
    our_decision[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  
  little_power[ind] = mean(little_decision, na.rm = TRUE)
  our_power[ind] = mean(our_decision, na.rm = TRUE)
  ind = ind+1
}

plot(seq(0.03, 0.3, length.out=7), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("bottomright",
       legend = c("Little's power", "Our power"),
       col = c("green", "blue"),
       pch = c(18, 18))

######### MAR ##############################
alpha = 0.05
n = 1000
MC = 100
d = 5

# Select the copula
cp = claytonCopula(param = c(1), dim = d)
# Generate the multivariate distribution (in this case it is just bivariate) with chisqal and t marginals
P = mvdc(copula = cp, margins = c(rep("norm",d)),
         paramMargins = rep(list(1),d) )

data = rMvdc(n, P)

xxx = seq(0.03, 0.3, length.out=7)

little_power = numeric(length = 7)
our_power = numeric(length = 7)
ind = 1
for (p in xxx){
  little_decision = logical(length = MC)
  our_decision = logical(length = MC)
  our_decision_P = logical(length = MC)
  for (i in 1:MC){
    X = delete_MAR_rank(data, p, c(1,2,5), cols_ctrl = c(3,4))
    
    little_decision[i] = mcar_test(data.frame(X))$p.value < alpha
    our_decision[i] = MCAR_meancovTest(X, alpha = 0.05, B = 99, type = "np")
  }
  
  little_power[ind] = mean(little_decision, na.rm = TRUE)
  our_power[ind] = mean(our_decision, na.rm = TRUE)
  ind = ind+1
}

plot(seq(0.03, 0.3, length.out=7), little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type="b")
lines(seq(0.03, 0.3, length.out=7), our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("bottomright",
       legend = c("Little's power", "Our power"),
       col = c("green", "blue"),
       pch = c(18, 18))
