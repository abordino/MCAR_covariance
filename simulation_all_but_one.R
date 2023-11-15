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
d = 5

little_power = c()
little_power_cov = c()
our_power = c()
R = c()

for(kkk in 1:10){
  
  #### POPULATION LEVEL ######
  SigmaS=list() 
  for(j in 1:d){
    A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
    tmp_0 = t(A) %*% A
    SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
  }
  
  v = 1:d
  pattern = list()
  for (i in 1:d){
    pattern[[i]] = v[-i]
  }
  result = computeR(patterns = pattern, SigmaS = SigmaS)
  
  R = c(R, result$R)
  
  ###### SAMPLE LEVEL, REPEATING THE TEST M TIMES #######
  little_decisions = c()
  little_decisions_cov = c()
  # our_decisions = c()
  for (i in 1:M){
    
    #### generate dataset
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:d){
      X[(1+(i-1)*n):(i*n), pattern[[i]]] = mvrnorm(n, rep(0, d-1), SigmaS[[i]])
    }
    X = as.matrix(X)
    
    ### run little's test
    little_decisions = c(little_decisions, little_test(X, alpha))
    little_decisions_cov = c(little_decisions_cov, little_test(X, alpha, "cov"))
    
    ### run our tests
    # our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
  }
  
  little_power = c(little_power, mean(little_decisions))
  little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  # our_power = c(our_power, mean(our_decisions))
}

plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "")
points(R, little_power_cov, col="orange", pch=18)
# points(R, our_power, col="blue", pch=18)
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))


######### it seems not to work, since the rotated data have R close to 1, even though 
######### they should mimic the null distribution. This is due to the fact that, due to 
######### sampling error, the sequence of correlation matrices is not compatible, but we
######### assume it's close to zero. For the cycle, and the block 3-cycle, this seems to be true,
######### while here, since we have many overlaps, R immediately takes big values! Even bigger than 
######### the original one.

#######################################################################################
################## R of at the pop level
SigmaS=list() 
for(j in 1:d){
  A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
  tmp_0 = t(A) %*% A
  SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
}

v = 1:d
pattern = list()
for (i in 1:d){
  pattern[[i]] = v[-i]
}
result = computeR(patterns = pattern, SigmaS = SigmaS)

R = result$R
R

##### generate data
X = data.frame(matrix(nrow = d*n, ncol = d))
for (i in 1:d){
  X[(1+(i-1)*n):(i*n), pattern[[i]]] = mvrnorm(n, rep(0, d-1), SigmaS[[i]])
}
X = as.matrix(X)


result = get_SigmaS(X)
SigmaS = result$SigmaS # this is what we called SigmaS_hat
patterns = result$pattern
n_pattern = result$n_pattern
data_pattern = result$data_pattern

for (i in length(SigmaS)){
  if (min(eigen(SigmaS[[i]])$values) < 10^-7){
    print("SigmaS is singular!")
    return(SigmaS)
  }
}

tmp = computeR(patterns, SigmaS)

R_hat_0 = tmp$R

Q_hat = tmp$Sigma/(1-tmp$R)
QS_hat = list()
for (i in 1:n_pattern){
  QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
}

#### rotate X, to make it look like it's from H0
rot_data_pattern = list()
for (i in 1:n_pattern){
  rot_data_pattern[[i]] = t((sqrtm(QS_hat[[i]])$B)%*%(solve(sqrtm(SigmaS[[i]])$B))%*%t(data_pattern[[i]]))
}

SS = list()
for (i in 1:n_pattern){
  SS[[i]] = cor(rot_data_pattern[[i]])
}

### observe that they don't differ to much, but it's enough to make it very incompatible
SS
QS_hat
