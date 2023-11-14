library(Rcsdp)
library(Matrix)

############## d-dimensional cycle #############
d=5
guess = 100
R = 0
while(abs(R - guess/2) > 0.00001){
  # R > 0
  SigmaS=list() 
  for(j in 1:d){
    A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
    tmp_0 = t(A) %*% A
    SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
  }
  
  # # R = 0
  # SigmaS=list() 
  # v = 1:d
  # A = matrix(runif(d^2)*2-1, ncol=d)
  # Sigma = cov2cor(t(A) %*% A)
  # for(j in 1:d){
  #   SigmaS[[j]]=as.matrix(Sigma[v[-j],v[-j]])
  # }
  
  
  # Matrix giving objective function (1/d)*tr(Sigma)
  C = c(list((1/d)*diag(d)), rep(list(matrix(0,d-1,d-1)),d))
  
  m=(d-1)*(d^2)/2 + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
  b=rep(0,m)
  A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,d-1,d-1)),d))
  
  for (k in 1:d){
    for (i in 1:(d-1)){
      for (j in i:(d-1)){
        b[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i] = SigmaS[[k]][i,j]
        
        e_i = rep(0, d-1); e_i[i]=1
        e_j = rep(0, d-1); e_j[j]=1
        A[[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i]] = Ablank
        A[[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i]][[k+1]] = e_i%*%t(e_j)/2 + e_j%*%t(e_i)/2
        
        v = 1:d
        tmp_1 = e_i%*%t(e_j)/2 +  e_j%*%t(e_i)/2
        tmp_2 = matrix(0,d,d)
        tmp_2[v[-k],v[-k]] = tmp_1
        A[[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i]][[1]] = tmp_2
      }
    }
  }
  
  for (i in 1:(d-1)){
    A[[(d-1)*(d^2)/2 + i]] = Ablank
    A[[(d-1)*(d^2)/2 + i]][[1]][1,1]=1
    A[[(d-1)*(d^2)/2 + i]][[1]][1+i,1+i]=-1
  }
  
  K = list(type=rep("s",d+1),size=c(d,rep(d-1,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  SDP=csdp(C, A, b, K) #Running the function
  R = 1-SDP$pobj
  
  ############### THE GUESS
  v = 1:d
  
  SigmaP=list() 
  for (j in 1:d){
    tmp_2 = matrix(0,d,d)
    tmp_2[v[-j],v[-j]] = SigmaS[[j]]
    SigmaP[[j]] = tmp_2
    
  }
  
  guess = 0
  best_i = 0; best_j = 0; best_k = 0; best_h = 0
  B = matrix(0,d,d)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      for (h in v[-c(i,j)]){
        for (k in v[-c(i,j,h)]){
          tmp = abs(SigmaP[[h]][i,j] - SigmaP[[k]][i,j])
          if (tmp > guess){
            guess = tmp
            best_i = i; best_j = j;best_k = k; best_h = h
          }
          if (tmp > B[i,j]){
            B[i,j] = tmp
          }
        }
      }
    }
  }
  
  matrix_means = matrix(0,d,d)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      for (h in v[-c(i,j)]){
          matrix_means[i,j] = matrix_means[i,j] + SigmaP[[h]][i,j]/(d-2)
        }
      }
    }

}

# THE RESULTS  
plot(0:1, 0:1, ylim=c(0,1), cex=0.1, main = paste("d = 3", " - simulation number ", ttt))
abline(h = R, col="red")
abline(h = guess, col="green")
best_i;best_j;best_h;best_k


# OPTIMAL SIGMA (dual)
SDP$X[[1]]


# OPTIMAL XS (primal)
YS=rep(list(matrix(0,d-1,d-1)),d)
ind = 1
for (k in 1:d){
  for (i in 1:(d-1)){
    for (j in i:(d-1)){
      tmp = as.numeric(abs(SDP$y[ind]/2) > 0.0001)*SDP$y[ind]
      YS[[k]][i,j] = tmp/2
      YS[[k]][j,i] = YS[[k]][j,i] + tmp/2
      ind = ind + 1
    }
  }
}
#
XS0 =rep(list(diag(1/(d-1),d-1)),d)
XS = c(list())
for (k in 1:d){
  XS[[k]] = d*YS[[k]] - XS0[[k]]
}


sum = 0
for (k in 1:d){
  sum = sum + -sum(diag(XS[[k]]%*%SigmaS[[k]]))/d
}

AstarXS = matrix(0,d,d)
XS_P=list()
for (j in 1:d){
  tmp_2 = matrix(0,d,d)
  tmp_2[v[-j],v[-j]] = XS[[j]]
  XS_P[[j]] = tmp_2
  AstarXS = AstarXS + XS_P[[j]]
}

XS



############## SIMULATION #############

for (d in 5:5){
  for (kkk in 1:1){
  # R > 0
  # SigmaS=list() 
  # for(j in 1:d){
  #   A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
  #   tmp_0 = t(A) %*% A
  #   SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
  # }
  
  # # R = 0
  # SigmaS=list()
  # v = 1:d
  # A = matrix(runif(d^2)*2-1, ncol=d)
  # Sigma = cov2cor(t(A) %*% A)
  # for(j in 1:d){
  #   SigmaS[[j]]=as.matrix(Sigma[v[-j],v[-j]])
  # }
    
    SigmaS=list()
    v = 1:d
    Sigma = diag(d)
    for(j in 1:d){
      SigmaS[[j]]=as.matrix(Sigma[v[-j],v[-j]])
    }
    
    for (i in 1:(d-2)){
      SigmaS[[d]][i,i+1] = runif(1, min = -1, max=1)/2
      SigmaS[[d]][i+1,i] = SigmaS[[1]][i,i+1]
    }
    
  
    
  
  # Matrix giving objective function (1/d)*tr(Sigma)
  C = c(list((1/d)*diag(d)), rep(list(matrix(0,d-1,d-1)),d))
  
  m=(d-1)*(d^2)/2 + d-1 # Number of constraints, all of form A(Sigma) + Z_\mathbb{S} = Sigma_\mathbb{S}
  b=rep(0,m)
  A=rep(list(0),m); Ablank=c(list(matrix(0,d,d)), rep(list(matrix(0,d-1,d-1)),d))
  
  for (k in 1:d){
    for (i in 1:(d-1)){
      for (j in i:(d-1)){
        b[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i] = SigmaS[[k]][i,j]
        
        e_i = rep(0, d-1); e_i[i]=1
        e_j = rep(0, d-1); e_j[j]=1
        A[[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i]] = Ablank
        A[[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i]][[k+1]] = e_i%*%t(e_j)/2 + e_j%*%t(e_i)/2
        
        v = 1:d
        tmp_1 = e_i%*%t(e_j)/2 +  e_j%*%t(e_i)/2
        tmp_2 = matrix(0,d,d)
        tmp_2[v[-k],v[-k]] = tmp_1
        A[[(k-1)*d*(d-1)/2 + (i-1)*d-(i-1)*i/2 + j+1-i]][[1]] = tmp_2
      }
    }
  }
  
  for (i in 1:(d-1)){
    A[[(d-1)*(d^2)/2 + i]] = Ablank
    A[[(d-1)*(d^2)/2 + i]][[1]][1,1]=1
    A[[(d-1)*(d^2)/2 + i]][[1]][1+i,1+i]=-1
  }
  
  K = list(type=rep("s",d+1),size=c(d,rep(d-1,d))) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  SDP=csdp(C, A, b, K) #Running the function
  R = 1-SDP$pobj
  
  ############### THE GUESS
  v = 1:d
  
  SigmaP=list() 
  for (j in 1:d){
    tmp_2 = matrix(0,d,d)
    tmp_2[v[-j],v[-j]] = SigmaS[[j]]
    SigmaP[[j]] = tmp_2
    
  }
  
  guess = 0
  best_i = 0; best_j = 0; best_k = 0; best_h = 0
  B = matrix(0,d,d)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      for (h in v[-c(i,j)]){
        for (k in v[-c(i,j,h)]){
          tmp = abs(SigmaP[[h]][i,j] - SigmaP[[k]][i,j])
          if (tmp > guess){
            guess = tmp
            best_i = i; best_j = j;best_k = k; best_h = h
          }
          if (tmp > B[i,j]){
            B[i,j] = tmp
          }
        }
      }
    }
  }
 
  plot(0:1, 0:1, ylim=c(0,2), cex=0.1, main = paste("d = ", d))
  abline(h = R, col="red")
  abline(h = guess, col="green")
  abline(h = guess/2, col = "blue")
  }
  
}


SDP$X[[1]]
SigmaP[[d]]
tilde_Sigma = SigmaP[[d]]/2 + diag((1-R)-1/2,d)
tilde_Sigma[1,1] = 1-R
eigen(tilde_Sigma)$values

