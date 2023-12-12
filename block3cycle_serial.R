source("computeR.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)

######### 3-cycle: setting 1 ############
alpha = 0.05
n = 500
MC = 1
d = 4

######## USING FOREACH AND DORNG

start.time = Sys.time()

# definition of P
P = 2*matrix(runif((d)^2)*2-1, ncol=d)/d

while(norm(P, type="2") > 0.5){
  P = matrix(runif((d)^2)*2-1, ncol=d)/d
}

normP = norm(P, type="2")
xxx = seq(1-2*normP^2, 0.9, length.out = 8)

little_power = numeric(length = 8)
little_power_cov = numeric(length = 8)
our_power = numeric(length = 8)
R = numeric(length = 8)

ind = 1
for (beta in xxx){
  print(beta)
  
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
  
  R[ind] = computeR(list(c(1:(2*d)),c(1:d,(2*d+1):(3*d)),c((d+1):(3*d))), SigmaS)$R

  little_decisions = logical(length = MC)
  little_decisions_cov = logical(length = MC)
  our_decisions = logical(length = MC)
  
  for (i in 1:MC){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X1 = mvrnorm(n, rep(0, 2*d), SigmaS[[1]])
    X2 = mvrnorm(n, rep(0, 2*d), SigmaS[[2]])
    X3 = mvrnorm(n, rep(0, 2*d), SigmaS[[3]])
    
    X = data.frame(matrix(nrow = 3*n, ncol = 3*d))
    X[1:n, c(1:(2*d))] = X1
    X[(n+1):(2*n), c(1:d,(2*d+1):(3*d))] = X2
    X[(2*n+1):(3*n), c((d+1):(3*d))] = X3
    X = as.matrix(X)
    
    ### run our tests
    little_decisions[i] = little_test(X, alpha = 0.05)
    little_decisions_cov[i] = little_test(X, alpha = 0.05, type = "cov")
    our_decisions[i] = MCAR_corr_test(X, alpha = 0.05, B = 99, type = "p")
  }
  
  little_power[ind] = mean(little_decisions)
  little_decisions_cov[ind] = mean(little_decisions_cov)
  our_power[ind] = mean(our_decisions)
  ind = ind+1
}

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

### plot the simulation

png(paste("block_3_cycle-",d, ".png"))
plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type = "b")
lines(R, little_power_cov, col="orange", pch=18, type = "b")
lines(R, our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("bottomright",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))
dev.off()