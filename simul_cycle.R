library(Rcsdp)
library(Matrix)

############## d-dimensional cycle #############
d=3
for (delta in seq(0.005, pi/2, length.out=50)){
R_rho =c()
XS_rho = list()
kkk = 1
for (rho in seq(0.01,pi/2,length.out=1000)){
  
if(pi - delta <= delta + rho){
  break
}
  
SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(pi - delta)
SigmaS[[1]][2,1] = cos(pi - delta)
SigmaS[[2]][1,2] = cos(delta)
SigmaS[[2]][2,1] = cos(delta)
SigmaS[[3]][1,2] = cos(rho)
SigmaS[[3]][2,1] = cos(rho)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
                                A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
R = 1-SDP$pobj

R_rho = c(R_rho,R)
############################
Sigma_opt = SDP$X[[1]] # Choice of Sigma with maximal trace satisfying SigmaS - A(Sigma) \succeq 0
qr(SDP$X[[1]])$rank # The rank of the upper block is
eigen(SDP$X[[1]])


############## Finding optimal choice of X_\mathbb{S} for what we call primal problem 
XS=rep(list(matrix(0,2,2)),d)
XS[[1]][1,1]=SDP$y[1]; XS[[d]][2,2]=SDP$y[2]; XS[[d]][1,2]=SDP$y[3]/2; XS[[d]][2,1]=SDP$y[3]/2
ind=3
for(j in 2:d){
  ind=ind+1; XS[[j]][1,1]=SDP$y[ind]
  ind=ind+1; XS[[j-1]][2,2]=SDP$y[ind]
  ind=ind+1; XS[[j-1]][1,2]=SDP$y[ind]/2; XS[[j-1]][2,1]=SDP$y[ind]/2
}

XS_rho[[kkk]] = XS
print(XS)
kkk = kkk + 1
}

plot(seq(0,pi/2,length.out=1000), R_rho, type = "l",col="red")
abline(h=(cos(delta)-cos(pi - delta))/2)
curve((cos(delta)-cos(pi - delta) - x)/2, add = T,col="blue")
curve((pi - 2*delta - x)*(cos(delta)-cos(pi - delta))/(2*pi-4*delta), add = T)
}

XS_new = list()
for (i in 1:3){
  XS_new[[i]] = 3*XS[[i]] - 0.5*diag(2)
}

delta=0.2
rho = 0
SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(pi - delta)
SigmaS[[1]][2,1] = cos(pi - delta)
SigmaS[[2]][1,2] = cos(delta)
SigmaS[[2]][2,1] = cos(delta)
SigmaS[[3]][1,2] = cos(rho)
SigmaS[[3]][2,1] = cos(rho)
# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
R = 1-SDP$pobj

############################
Sigma_opt = SDP$X[[1]] # Choice of Sigma with maximal trace satisfying SigmaS - A(Sigma) \succeq 0
qr(SDP$X[[1]])$rank # The rank of the upper block is
eigen(SDP$X[[1]])


############## Finding optimal choice of X_\mathbb{S} for what we call primal problem 
XS=rep(list(matrix(0,2,2)),d)
XS[[1]][1,1]=SDP$y[1]; XS[[d]][2,2]=SDP$y[2]; XS[[d]][1,2]=SDP$y[3]/2; XS[[d]][2,1]=SDP$y[3]/2
ind=3
for(j in 2:d){
  ind=ind+1; XS[[j]][1,1]=SDP$y[ind]
  ind=ind+1; XS[[j-1]][2,2]=SDP$y[ind]
  ind=ind+1; XS[[j-1]][1,2]=SDP$y[ind]/2; XS[[j-1]][2,1]=SDP$y[ind]/2
}

for (i in 1:3){
  XS[[i]] = 3*XS[[i]] - 0.5*diag(2)
}

AstarXS=matrix(0,d,d)
AstarXS[1,1]=XS[[1]][1,1]+XS[[d]][2,2]
AstarXS[1,d]=XS[[d]][1,2]; AstarXS[d,1]=AstarXS[1,d]
for(j in 2:d){
  AstarXS[j,j]=XS[[j-1]][2,2]+XS[[j]][1,1]
  AstarXS[j,j-1]=XS[[j-1]][1,2]; AstarXS[j-1,j]=AstarXS[j,j-1]
}

#################### \mathbb{S} = {{1,2},{1,3}} ######################

NSim=20000
Rind=rep(0,NSim)
incons=rep(0,NSim)
for(nsim in 1:NSim){
  # Input Sigma_\mathbb{S} should be a list with elements Sigma_{1,2},Sigma_{1,3}
  SigmaS=list() #Random 2x2 PSD matrices with trace standardisation
  for(j in 1:2){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=x%*%t(x) + y%*%t(y)
  }
  barTrace=(1/2)*(SigmaS[[1]][1,1]+SigmaS[[2]][1,1])+SigmaS[[1]][2,2]+SigmaS[[2]][2,2]
  SigmaS[[1]]=SigmaS[[1]]*3/barTrace; SigmaS[[2]]=SigmaS[[2]]*3/barTrace
  
  # Matrix giving objective function (1/3)*tr(Sigma)
  C = c(list((1/3)*diag(3)), rep(list(matrix(0,2,2)),2) )
  
  m=6  # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
  b=rep(0,m)
  A=rep(list(0),m); Ablank=c(list(matrix(0,3,3)), rep(list(matrix(0,2,2)),2))
  b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
  b[2]=SigmaS[[2]][1,1]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[3]][1,1]=1 # Constraint for {1,3},1,1
  b[3]=SigmaS[[1]][2,2]; A[[3]]=Ablank; A[[3]][[1]][2,2]=1; A[[3]][[2]][2,2]=1 # Constraint for {1,2},2,2
  b[4]=SigmaS[[2]][2,2]; A[[4]]=Ablank; A[[4]][[1]][3,3]=1; A[[4]][[3]][2,2]=1 # Constraint for {1,3},3,3
  b[5]=SigmaS[[1]][1,2]; A[[5]]=Ablank; A[[5]][[1]][1,2]=1/2; A[[5]][[1]][2,1]=1/2 
  A[[5]][[2]][1,2]=1/2; A[[5]][[2]][2,1]=1/2 # Constraint for {1,2},1,2
  b[6]=SigmaS[[2]][1,2]; A[[6]]=Ablank; A[[6]][[1]][1,3]=1/2; A[[6]][[1]][3,1]=1/2 
  A[[6]][[3]][1,2]=1/2; A[[6]][[3]][2,1]=1/2 # Constraint for {1,3},1,3
  
  K = list(type=rep("s",3),size=c(3,rep(2,2))) # Variables are the PSD matrices Sigma,Z_{1,2},Z_{1,3}
  
  SDP=csdp(C, A, b, K) #Running the function
  
  SDP$X[[1]] # Choice of Sigma with maximal trace satisfying SigmaS - A(Sigma) \succeq 0
  
  Rind[nsim]=1-SDP$pobj # Value of R(P_\mathbb{S}) (computed using what we call dual form in paper)
  incons[nsim]=abs(SigmaS[[1]][1,1]-SigmaS[[2]][1,1])
}
plot(Rind,incons)

#################### same for 4-cycle ######################
NSim=20000
Rind=rep(0,NSim)
incons=rep(0,NSim)
for(nsim in 1:NSim){
  d=4
  
  # Input Sigma_\mathbb{S} should be a list with elements Sigma_{1,2},...,Sigma_{d-1,d},Sigma_{d,1}
  SigmaS=list() #Random 2x2 PSD matrices
  for(j in 1:d){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=x%*%t(x) + y%*%t(y)
  }
  
  SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
  for(j in 1:d){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }
  
  # # user defined matrices
  # rho = runif(1,0,pi)
  # SigmaS[[1]][1,2] = -1
  # SigmaS[[1]][2,1] = -1
  # SigmaS[[2]][1,2] = -1
  # SigmaS[[2]][2,1] = -1
  # SigmaS[[3]][1,2] = -1
  # SigmaS[[3]][2,1] = -1
  # SigmaS[[4]][1,2] = cos(rho)
  # SigmaS[[4]][2,1] = cos(rho)
  
  # Matrix giving objective function (1/d)*tr(Sigma)
  C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))
  
  m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
  b=rep(0,m)
  A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
  b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
  b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
  b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
  A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
  ind=3
  for(j1 in 2:d){
    ind=ind+1
    b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
    ind=ind+1
    b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
    ind=ind+1
    b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
    A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
    
  }
  
  for (i in 1:(d-1)){
    A[[3*d + i]] = Ablank
    A[[3*d + i]][[1]][1,1]=1
    A[[3*d + i]][[1]][1+i,1+i]=-1
  }
  
  K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  SDP=csdp(C, A, b, K) #Running the function
  
  ############## Check if KKT conditions hold true
  a = SigmaS[[1]][1,2];b = SigmaS[[2]][1,2];c = SigmaS[[3]][1,2];d = SigmaS[[4]][1,2]
  specified = sort(c(a,b,c,d))
  theta = acos(specified)
  
  Rind[nsim]=1-SDP$pobj #- sin(max(0, (theta[1]+theta[2]+theta[3]-theta[4]-2*pi)/2) + (max(0, (theta[1]-(theta[2]+theta[3]+theta[4]))))/2)
  # Value of R(P_\mathbb{S}) (computed using what we call dual form in paper)
  
  # incons[nsim]= max(0, (theta[1]+theta[2]+theta[3]-theta[4]-2*pi)) + 
  #   max(0, (theta[1]-(theta[2]+theta[3]+theta[4])))
  
  
  ### general case
  # num = (b*c - a*d)*(b*d-a*c)*(-a*b+c*d)
  # den1 = (-a-b+c+d)*(a-b+c-d)*(a-b-c+d)*(a+b+c+d)
  # if (max(0, (theta[1]+theta[2]+theta[3]-theta[4]-2*pi)) + 
  #     max(0, (theta[1]-(theta[2]+theta[3]+theta[4]))) > 0){
  #   incons[nsim] = max(0,-1+2*num/sqrt(den1*num))/(2*sqrt(2)-1)
  # }
  # else{
  #   incons[nsim] = 0
  # }
  
  
  # incons[nsim]= rho
}

plot(incons, Rind, pch = 20, cex = 0.1)
curve(cos(x/2)^2, col = "green", add = T)
curve(0*x, col = "red", add = T)
curve(sin(x/2)^(9/10), col = "red", add = T)
curve(sin(x/2)^(1/2), col = "green", add = T)
curve(sin(x/2)^(2), col = "blue", add = T)
curve(cos(x/2)^(2), col = "blue", add = T)






########## REDUCTION 3 ###########
library(Rcsdp)
library(Matrix)
d = 3
R_rho = c()
  for (rho in seq(from=0, to=pi, length.out=1000)){
    # Input Sigma_\mathbb{S} should be a list with elements Sigma_{1,2},...,Sigma_{d-1,d},Sigma_{d,1}
    
    SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
    for(j in 1:d){
      x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
    }
    
    for (i in 1:(d-1)){
      SigmaS[[i]][1,2] = 0
      SigmaS[[i]][2,1] = 0
    }
    
    SigmaS[[d]][1,2] = cos(rho)
    SigmaS[[d]][2,1] = cos(rho)
    
    
    # Matrix giving objective function (1/d)*tr(Sigma)
    C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))
    
    m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
    b=rep(0,m)
    A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
    b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
    b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
    b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
    A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
    ind=3
    for(j1 in 2:d){
      ind=ind+1
      b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
      ind=ind+1
      b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
      ind=ind+1
      b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
      A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
      
    }
    
    for (i in 1:(d-1)){
      A[[3*d + i]] = Ablank
      A[[3*d + i]][[1]][1,1]=1
      A[[3*d + i]][[1]][1+i,1+i]=-1
    }
    
    K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
    
    SDP=csdp(C, A, b, K) #Running the function
    R = 1-SDP$pobj
    R_rho = c(R_rho, R)
    
  }

plot(seq(from=0, to=pi, length.out=1000), R_rho, ylim=c(0,1))

#########################################################
####### New conjecture
library(Rcsdp)
library(Matrix)

d=3

for (kkk in 1:100000){
SigmaS=list() 
for(j in 1:3){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

theta_2 = runif(1,0,pi/2)
theta_3 = runif(1,0,pi/2)
theta_1 = runif(1,theta_2+theta_3,pi)

SigmaS[[1]][1,2] = cos(theta_1)
SigmaS[[1]][2,1] = cos(theta_1)
SigmaS[[2]][1,2] = cos(theta_2)
SigmaS[[2]][2,1] = cos(theta_2)
SigmaS[[3]][1,2] = cos(theta_3)
SigmaS[[3]][2,1] = cos(theta_3)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
R = 1-SDP$pobj

h = (theta_1-theta_2-theta_3)/(theta_1-theta_2)
Conj = h*(cos(theta_2)-cos(theta_1))/2

# plot(c(0,1),c(0,1))
# abline(h=R, col='red')
# abline(h=Conj, col='blue')

  if (Conj > R){
    print("not true!!!!!")
    break
  }
}

############## Finding optimal choice of X_\mathbb{S} for what we call primal problem 
XS=rep(list(matrix(0,2,2)),d)
XS[[1]][1,1]=SDP$y[1]; XS[[d]][2,2]=SDP$y[2]; XS[[d]][1,2]=SDP$y[3]/2; XS[[d]][2,1]=SDP$y[3]/2
ind=3
for(j in 2:d){
  ind=ind+1; XS[[j]][1,1]=SDP$y[ind]
  ind=ind+1; XS[[j-1]][2,2]=SDP$y[ind]
  ind=ind+1; XS[[j-1]][1,2]=SDP$y[ind]/2; XS[[j-1]][2,1]=SDP$y[ind]/2
}
XS
h
XS_new = list()
for (i in 1:3){
  XS_new[[i]] = 3*XS[[i]] - 0.5*diag(2)
}

SigmaS
1 - sum(diag(XS[[1]]%*%SigmaS[[1]])) - sum(diag(XS[[2]]%*%SigmaS[[2]])) - sum(diag(XS[[3]]%*%SigmaS[[3]]))

stop = F
for (c_1 in seq(0.254671,0.3,length.out=10000)){
  for (c_2 in seq(3.0997,3.1,length.out=10000)){
    YS = list()
    YS[[1]] = matrix(c(c_1,c_1,c_1,c_1), nrow = 2)
    YS[[2]] = matrix(c(c_1,-c_1,-c_1,c_1), nrow = 2)
    YS[[3]] = matrix(c(c_2,-c_2,-c_2,c_2), nrow = 2)
    trial = 1 - sum(diag(YS[[1]]%*%SigmaS[[1]])) - sum(diag(YS[[2]]%*%SigmaS[[2]])) - sum(diag(YS[[3]]%*%SigmaS[[3]]))
    print(c_1)
    print(c_2)
    print(trial)
    
    if(abs(trial-Conj)<0.00001 && trial>0){
      print("maybe found")
      stop = T
      break
    }
  }
  if (stop){
    break
  }
}

stop = F
  for (c_2 in seq(148.8045,200,length.out=10000000)){
    YS = list()
    YS[[1]] = matrix(c(h/4,h/4,h/4,h/4), nrow = 2)
    YS[[2]] = matrix(c(h/4,-h/4,-h/4,h/4), nrow = 2)
    YS[[3]] = matrix(c(c_2,-c_2,-c_2,c_2), nrow = 2)
    trial = 1 - sum(diag(YS[[1]]%*%SigmaS[[1]])) - sum(diag(YS[[2]]%*%SigmaS[[2]])) - sum(diag(YS[[3]]%*%SigmaS[[3]]))
    print(c_2)
    print(trial)
    
    if(abs(trial-Conj)<0.1 && trial>0){
      print("maybe found")
      stop = T
      break
    }
  }

XS = list()
for (i in 1:3){
  XS[[i]] = 3*YS[[i]] - 0.5*diag(2)
}

AstarXS=matrix(0,3,3)
AstarXS[1,1]=XS[[1]][1,1]+XS[[3]][2,2]
AstarXS[1,3]=XS[[3]][1,2]; AstarXS[3,1]=AstarXS[1,3]
for(j in 2:3){
  AstarXS[j,j]=XS[[j-1]][2,2]+XS[[j]][1,1]
  AstarXS[j,j-1]=XS[[j-1]][1,2]; AstarXS[j-1,j]=AstarXS[j,j-1]
}

Y = 0*diag(3)
for (e in seq(0,3,length.out=100000)){
  Y[1,1] = -e/2; Y[3,3] = -e/2; Y[2,2] = e
  v = eigen(AstarXS+Y)$values
  print("---------------------------")
  print(e)
  print(v)
  if(v[1] >= 0 && v[2] >= 0 && v[3] >= 0){
    print("Found Y!!!")
    break
  }
}




####### New conjecture
library(Rcsdp)
library(Matrix)

d=3

theta_2 = runif(1,0,pi/2)
theta_1 = runif(1,theta_2,pi)
while(theta_1 - theta_2 > pi/2){
  theta_2 = runif(1,0,pi/2)
  theta_1 = runif(1,theta_2,pi)
}

R_theta_3 = c()
cos_theta_3 = c()
s1 = c()
s2 = c()
s3 = c()
for (theta_3 in seq(0,(theta_1-theta_2),length.out=1000)){
  SigmaS=list() 
  for(j in 1:3){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }
  
  SigmaS[[1]][1,2] = cos(theta_1)
  SigmaS[[1]][2,1] = cos(theta_1)
  SigmaS[[2]][1,2] = cos(theta_2)
  SigmaS[[2]][2,1] = cos(theta_2)
  SigmaS[[3]][1,2] = cos(theta_3)
  SigmaS[[3]][2,1] = cos(theta_3)
  
  # Matrix giving objective function (1/d)*tr(Sigma)
  C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))
  
  m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
  b=rep(0,m)
  A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
  b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
  b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
  b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
  A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
  ind=3
  for(j1 in 2:d){
    ind=ind+1
    b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
    ind=ind+1
    b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
    ind=ind+1
    b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
    A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
    
  }
  
  for (i in 1:(d-1)){
    A[[3*d + i]] = Ablank
    A[[3*d + i]][[1]][1,1]=1
    A[[3*d + i]][[1]][1+i,1+i]=-1
  }
  
  K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  SDP=csdp(C, A, b, K) #Running the function
  R = 1-SDP$pobj
  R_theta_3 = c(R_theta_3, R)
  
  cos_phi_3 = SDP$X[[1]][1,3]/(1-R)
  cos_theta_3 = c(cos_theta_3, cos_phi_3)
  
  tmp = (SigmaS[[1]][1,2] - SDP$X[[1]][1,2])/R
  s1 = c(s1, tmp)
  
  tmp = (SigmaS[[2]][1,2] - SDP$X[[2]][1,2])/R
  s2 = c(s2, tmp)
  
  tmp = (SigmaS[[3]][1,2] - SDP$X[[3]][1,2])/R
  s3 = c(s3, tmp)
}

plot(seq(0,(theta_1-theta_2),length.out=1000), R_theta_3, col="red", type="l", 
     xlim=c(0,0.01),ylim = c(h=0.5*(cos(theta_2)-cos(theta_1))-0.0001,h=0.5*(cos(theta_2)-cos(theta_1))))
abline(h=0.5*(cos(theta_2)-cos(theta_1)))
curve(0.5*(cos(theta_2)-cos(theta_1))*(theta_1-theta_2-x)/(theta_1-theta_2),add=T,col="blue")

plot(seq(0,(theta_1-theta_2),length.out=1000), R_theta_3, col="red", type="l")
abline(h=0.5*(cos(theta_2)-cos(theta_1)))
curve(0.5*(cos(theta_2)-cos(theta_1))*(theta_1-theta_2-x)/(theta_1-theta_2),add=T,col="blue")

plot(seq(0,(theta_1-theta_2),length.out=1000), cos_theta_3, col="red", type="l")

plot(seq(0,(theta_1-theta_2),length.out=1000), s3, col="red", type="l")
plot(seq(0,(theta_1-theta_2),length.out=1000), s2, col="red", type="l")
plot(seq(0,(theta_1-theta_2),length.out=1000), s1, col="red", type="l")
theta_1
theta_2
theta_3

plot(seq(0,(theta_1-theta_2),length.out=1000), 1-cos_theta_3, col="red", type="l")
curve(1-cos(x),add=T, col="blue")

plot(seq(0,(theta_1-theta_2),length.out=1000), acos(cos_theta_3), col="red", type="l")
curve(x^1,add=T, col="blue")

plot(seq(0,(theta_1-theta_2),length.out=1000), 
     (1-cos(seq(0,(theta_1-theta_2),length.out=1000)))/(1-cos_theta_3), col="orange", type="l",
     xlim = c(0, theta_1-theta_2), ylim = c(1-R_theta_3[1],1))

plot(seq(0,(theta_1-theta_2),length.out=1000), 1-R_theta_3, col="red", type="l")


############# general d
d=6

for (kkk in 1:100000){
  SigmaS=list() 
  for(j in 1:d){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }
  
  theta = c(runif(1,0,pi))
  for (i in 2:d){
    theta = c(theta, runif(1,0,theta[i-1]))
  }
  
  for (i in 1:d){
    SigmaS[[i]][1,2] = cos(theta[i])
    SigmaS[[i]][2,1] = cos(theta[i])
  }
  
  # Matrix giving objective function (1/d)*tr(Sigma)
  C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))
  
  m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
  b=rep(0,m)
  A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
  b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
  b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
  b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
  A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
  ind=3
  for(j1 in 2:d){
    ind=ind+1
    b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
    ind=ind+1
    b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
    ind=ind+1
    b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
    A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
    
  }
  
  for (i in 1:(d-1)){
    A[[3*d + i]] = Ablank
    A[[3*d + i]][[1]][1,1]=1
    A[[3*d + i]][[1]][1+i,1+i]=-1
  }
  
  K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  SDP=csdp(C, A, b, K) #Running the function
  R = 1-SDP$pobj
  
  
  h = sum(theta*c(1, rep(-1,d-1)))/(theta[1]-theta[2])
  Conj = h*(cos(theta[2])-cos(theta[1]))/2
  
  # plot(c(0,1),c(0,1))
  # abline(h=R, col='red')
  # abline(h=Conj, col='blue')
  
  if (Conj > R){
    print("not true!!!!!")
    print(R)
    print(Conj)
    print(theta)
    break
  }
}



d=5

R_theta5 = c()
R_theta5_prime = c()
R_mean = c()

for (kkk in 1:1000){
  
theta_2 = runif(1,0,pi/2)
theta_1 = runif(1,theta_2,pi)
while(theta_1 - theta_2 > pi/2){
  theta_2 = runif(1,0,pi/2)
  theta_1 = runif(1,theta_2,pi)
}
theta_3 = runif(1,0, theta_1-theta_2)
theta_4 = runif(1,0, theta_1-theta_2-theta_3)
theta_5 = runif(1,0, theta_1-theta_2-theta_3-theta_4)

#################### R theta_5
SigmaS=list() 
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(theta_1)
SigmaS[[1]][2,1] = cos(theta_1)
SigmaS[[2]][1,2] = cos(theta_2)
SigmaS[[2]][2,1] = cos(theta_2)
SigmaS[[3]][1,2] = cos(theta_3)
SigmaS[[3]][2,1] = cos(theta_3)
SigmaS[[4]][1,2] = cos(theta_4)
SigmaS[[4]][2,1] = cos(theta_4)
SigmaS[[5]][1,2] = cos(theta_5)
SigmaS[[5]][2,1] = cos(theta_5)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
R = 1-SDP$pobj


#################### R theta_5_prime
theta_5_prime = runif(1,0, theta_1-theta_2-theta_3-theta_4)

SigmaS[[5]][1,2] = cos(theta_5)
SigmaS[[5]][2,1] = cos(theta_5)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
tmp = 1-SDP$pobj
R_theta5 = c(R_theta5, tmp)

#################### R theta_5
SigmaS=list() 
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(theta_1)
SigmaS[[1]][2,1] = cos(theta_1)
SigmaS[[2]][1,2] = cos(theta_2)
SigmaS[[2]][2,1] = cos(theta_2)
SigmaS[[3]][1,2] = cos(theta_3)
SigmaS[[3]][2,1] = cos(theta_3)
SigmaS[[4]][1,2] = cos(theta_4)
SigmaS[[4]][2,1] = cos(theta_4)
SigmaS[[5]][1,2] = cos(theta_5)
SigmaS[[5]][2,1] = cos(theta_5)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
tmp = 1-SDP$pobj
R_theta5_prime = c(R_theta5_prime, tmp)


#################### R means

SigmaS[[5]][1,2] = cos((theta_5+theta_5_prime))/2
SigmaS[[5]][2,1] = cos((theta_5+theta_5_prime))/2

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
tmp = 1-SDP$pobj
R_mean = c(R_mean, tmp)
}

plot(R_mean, R_theta5+R_theta5_prime, col="red")




###########################################
theta_2 = runif(1,0,pi/2)
theta_1 = runif(1,theta_2,pi)
while(theta_1 - theta_2 > pi/2){
  theta_2 = runif(1,0,pi/2)
  theta_1 = runif(1,theta_2,pi)
}
theta_3 = runif(1,0, theta_1-theta_2)
theta_4 = runif(1,0, theta_1-theta_2-theta_3)
theta_5 = runif(1,0, theta_1-theta_2-theta_3-theta_4)

#################### R_Sigma_5
d=5
SigmaS=list() 
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(theta_1)
SigmaS[[1]][2,1] = cos(theta_1)
SigmaS[[2]][1,2] = cos(theta_2)
SigmaS[[2]][2,1] = cos(theta_2)
SigmaS[[3]][1,2] = cos(theta_3)
SigmaS[[3]][2,1] = cos(theta_3)
SigmaS[[4]][1,2] = cos(theta_4)
SigmaS[[4]][2,1] = cos(theta_4)
SigmaS[[5]][1,2] = cos(theta_5)
SigmaS[[5]][2,1] = cos(theta_5)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
R = 1-SDP$pobj
Sigma = SDP$X[[1]]


RE = c()
RB = c()
for(phi in seq(theta_4+theta_5, UB, length.out=1000)){
  
####### R_E
d = 3
SigmaS=list() 
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(theta_4)
SigmaS[[1]][2,1] = cos(theta_4)
SigmaS[[2]][1,2] = cos(theta_5)
SigmaS[[2]][2,1] = cos(theta_5)
SigmaS[[3]][1,2] = cos(phi)
SigmaS[[3]][2,1] = cos(phi)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
tmp = 1-SDP$pobj
RE = c(RE, tmp)

####### R_B
d = 4
SigmaS=list() 
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(phi)
SigmaS[[1]][2,1] = cos(phi)
SigmaS[[2]][1,2] = cos(theta_1)
SigmaS[[2]][2,1] = cos(theta_1)
SigmaS[[3]][1,2] = cos(theta_2)
SigmaS[[3]][2,1] = cos(theta_2)
SigmaS[[4]][1,2] = cos(theta_3)
SigmaS[[4]][2,1] = cos(theta_3)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
tmp = 1-SDP$pobj
RB = c(RB, tmp)
}

plot(seq(theta_4+theta_5, UB, length.out=1000), RE, type = "l", col="pink",
     xlim=c(theta_4+theta_5, UB), ylim=c(0,R+0.0001))
lines(seq(theta_4+theta_5, UB, length.out=1000), RB, type = "l", col="orange")
lines(seq(theta_4+theta_5, UB, length.out=1000), RE+RB, type = "l", col="green")
abline(h=R, col="red")



###########################################
library(Rcsdp)
theta_2 = runif(1,0,pi/2)
theta_1 = runif(1,theta_2,pi)
while(theta_1 - theta_2 > pi/2){
  theta_2 = runif(1,0,pi/2)
  theta_1 = runif(1,theta_2,pi)
}
theta = c(theta_1, theta_2)
theta_3 = runif(1,0, theta_1-theta_2)
theta = c(theta, tmp)
theta_4 = runif(1,0, theta_1-theta_2-theta_3)
theta = c(theta, tmp)
theta_5 = runif(1,0, theta_1-theta_2-theta_3-theta_4)
theta = c(theta, tmp)

sum(theta*(c(1, rep(-1,4))))


#################### R_Sigma_5
d=5
SigmaS=list() 
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(theta_1)
SigmaS[[1]][2,1] = cos(theta_1)
SigmaS[[2]][1,2] = cos(theta_2)
SigmaS[[2]][2,1] = cos(theta_2)
SigmaS[[3]][1,2] = cos(theta_3)
SigmaS[[3]][2,1] = cos(theta_3)
SigmaS[[4]][1,2] = cos(theta_4)
SigmaS[[4]][2,1] = cos(theta_4)
SigmaS[[5]][1,2] = cos(theta_5)
SigmaS[[5]][2,1] = cos(theta_5)

# Matrix giving objective function (1/d)*tr(Sigma)
C = c(list((1/d)*diag(d)), rep(list(matrix(0,2,2)),d))

m=3*d + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
b=rep(0,m)
A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,2,2)),d))
b[1]=SigmaS[[1]][1,1]; A[[1]]=Ablank; A[[1]][[1]][1,1]=1; A[[1]][[2]][1,1]=1 # Constraint for {1,2},1,1
b[2]=SigmaS[[d]][2,2]; A[[2]]=Ablank; A[[2]][[1]][1,1]=1; A[[2]][[d+1]][2,2]=1 # Constraint for {d,1},1,1
b[3]=SigmaS[[d]][1,2]; A[[3]]=Ablank; A[[3]][[1]][1,d]=1/2; A[[3]][[1]][d,1]=1/2 
A[[3]][[d+1]][1,2]=1/2; A[[3]][[d+1]][2,1]=1/2 # Constraint for {d,1},d,1
ind=3
for(j1 in 2:d){
  ind=ind+1
  b[ind]=SigmaS[[j1]][1,1]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1+1]][1,1]=1 # Constraint for {j1,j1+1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][2,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1]=1; A[[ind]][[j1]][2,2]=1 # Constraint for {j1-1,j1},j1,j1
  ind=ind+1
  b[ind]=SigmaS[[j1-1]][1,2]; A[[ind]]=Ablank; A[[ind]][[1]][j1,j1-1]=1/2; A[[ind]][[1]][j1-1,j1]=1/2 
  A[[ind]][[j1]][1,2]=1/2; A[[ind]][[j1]][2,1]=1/2 # Constraint for {j1-1,j1},j1-1,j1
  
}

for (i in 1:(d-1)){
  A[[3*d + i]] = Ablank
  A[[3*d + i]][[1]][1,1]=1
  A[[3*d + i]][[1]][1+i,1+i]=-1
}

K = list(type=rep("s",d+1),size=c(d,rep(2,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}

SDP=csdp(C, A, b, K) #Running the function
R = 1-SDP$pobj
Sigma = SDP$X[[1]]/(1-R)

phi = c(acos(Sigma[1,2]),acos(Sigma[3,2]),acos(Sigma[3,4]),acos(Sigma[4,5]),
        acos(Sigma[1,5]))
sum(phi*(c(1, rep(-1,4))))

