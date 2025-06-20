rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Desktop/simulation_regularised/")

# Load external functions and libraries
source("MCARtest/computeR.R")
source("MCARtest/find_SigmaS.R")
source("MCARtest/indexConsistency.R")
source("bootstrap_test_cycle.R")
source("little_test.R")

library(missMethods)
library(norm)
library(naniar)
library(MASS)
library(latex2exp)
library(compositions)

# #----------------------------------------------------------------------------------------
# ######### d-cycle: high-dimensional ############
# #----------------------------------------------------------------------------------------
set.seed(180966)

d = 200
n = 200
MC = 500

alpha = 0.05
angle = pi/(2*(d-1))

our_power_corr = c()
our_power = c()
R = c()

for(t1 in seq(pi/2, 5*pi/8, length.out = 8)){
  print("------------------------------------------------------------------------")
  print(t1)
  print("------------------------------------------------------------------------")
  
  #----------------------------------------------------------------------------------------
  #### POPULATION LEVEL ######
  #----------------------------------------------------------------------------------------
  SigmaS=list()
  for(j in 1:d){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
    SigmaS[[j]][1,2] = cos(angle)
    SigmaS[[j]][2,1] = cos(angle)
  }
  
  SigmaS[[1]][1,2] = cos(t1)
  SigmaS[[1]][2,1] = cos(t1)
  
  patterns = list()
  v = 1:d
  for (i in 1:(d-1)){
    patterns[[i]] = c(v[i], v[i+1])
  }
  patterns[[d]] = c(1,d)
  R = c(R, computeR(patterns, SigmaS)$R)
  
  
  #----------------------------------------------------------------------------------------
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  #----------------------------------------------------------------------------------------
  our_decisions = c()
  our_decisions_corr = c()
  for (i in 1:MC){
    
    #----------------------------------------------------------------------------------------
    #### generate dataset from patter S = {{1,2},{2,3}, ..., {1,d}}
    #----------------------------------------------------------------------------------------
    d0 = d; n0 = n*d; patterns0 = patterns
    X = data.frame(matrix(nrow = n0, ncol = d0))
    
    data_pattern0 = list()
    mu_S0 = list(); C_S0 = list(); sigma_squared_S0 = list(); SigmaS0 = list()
    for (i in 1:d){
      tmp = mvrnorm(n, rep(0, 2), SigmaS[[i]])
      X[(1+(i-1)*n):(i*n), patterns0[[i]]] = tmp
      data_pattern0[[i]] = tmp
      
      mu_S0[[i]] = colMeans(data_pattern0[[i]])
      C_S0[[i]] = var(data_pattern0[[i]])
      sigma_squared_S0[[i]] = diag(C_S0[[i]])
      SigmaS0[[i]] = cov2cor(C_S0[[i]])
    }
    
    X = as.matrix(X)
    
    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    p_R = corr.compTest.cycle(X, B=99)
    p_M = mean.consTest.cycle(X, B=99)
    p_V = var.consTest.cycle(X, B=99)
    
    our_decisions = c(our_decisions, max(c(p_R, p_M, p_V) < alpha/3))
    our_decisions_corr = c(our_decisions_corr, p_R < alpha)
    
  }
  
  our_power = c(our_power, mean(our_decisions))
  our_power_corr = c(our_power_corr, mean(our_decisions_corr))
}

#----------------------------------------------------------------------------------------
# ######### save data ############
#----------------------------------------------------------------------------------------
data_to_save = data.frame(R, our_power, our_power_corr, alpha)
write.csv(data_to_save, "simulCycle/results/200_cyclePar.csv", row.names = FALSE)

#----------------------------------------------------------------------------------------
# ######### plot results ############
#----------------------------------------------------------------------------------------

png("simulCycle/pictures/200_cyclePar.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(R, our_power, col="blue", ylim = c(0,1), pch=21,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylab = "Power", type = "b")
lines(R, our_power_corr, col="darkviolet", pch=25, type = "b")
lines(R, rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.4,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("Omnibus", TeX(r'($p_R$)')),
       col = c("blue", "darkviolet"),
       pch = c(21, 25))
dev.off()



