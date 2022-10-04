LOCUS <- function(Y, q, V, MaxIteration=100, penalty="SCAD", phi = 0.9,
                  approximation=TRUE, preprocess=TRUE, 
                  espli1=0.001, espli2=0.001, rho=0.95, silent=FALSE)
  # Y: connectivity data of dimension N x K, N is number of subjects, K is number of edges. 
  # q: Number of subnetworks to extract.
  # V: Number of nodes in network.
  # MaxIteration = 100: number of maximum iterations.
  # espli1=0.001, espli2=0.001: toleration for convergence on S and A. 
  # penalty = "SCAD": sparsity penalty, this can be NULL, SCAD, or L1.
  # phi = 0.9: tuning parameter for penalty
  # rho = 0.95: tuning parameter for selecting number of ranks in each subnetwork's decomposition. 
  # silent: whether to print intermediate results. 
  # preprocess: whether to preprocess the data Y (dont change if you are not sure) 
  # approximation = TRUE: whether to use approximated algorithm based on SVD (don't change if you are not sure).
{
  # demean the data
  Y = sweep(Y,2,apply(Y,2,mean),"-") 
  # preprocess the data if True
  if(preprocess)
  { 
    Yraw = Y
    Y = Locus_preprocess(Y,q) 
  }
  
  K = dim(Y)[2]          # Number of edges
  if(V != (sqrt(1+8*K)+1)/2)
  {
    stop("V is not correctly specified! Please double check the dimension of your input data.")
  }
  N = dim(Y)[1]          # Number of subjects (ICs)
  
  # Initial Estimation
  theta_ini = Locus_initial(Y,q,V,rho=rho)
  A = theta_ini$A; S = theta_ini$S
  theta = theta_ini$theta
  
  # Update Parameters
  Iter = 1
  while(Iter <= MaxIteration)
  {
    if(approximation)
    {
      theta_new = Locus_update_approx(Y,A,theta,penalt= penalty,lambda_ch = phi, gamma = 2.1,imput_method = "Previous",silent = silent)
    }else
    {
      theta_new = Locus_update(Y,A,theta,penalt= penalty,lambda_ch = phi, gamma = 2.1,silent = silent)
    }
    
    # orthogonize A here
    if(preprocess){
      A_new = far::orthonormalization(theta_new$A)
    }else{
      A_new = theta_new$A
    }
    
    S_new = theta_new$S; theta_new  = theta_new$theta
    
    errS = norm(as.matrix(S_new-S))/norm(as.matrix(S))
    errA = norm(as.matrix(A_new-A))/norm(as.matrix(A))
    
    if(sum(is.na(c(errS,errA)))>0){return(list(Conver=FALSE,A=A,S=S,theta=theta))}
    
    if(!silent)
    {
      message(paste("Iter ",Iter,"; Percentage change on S: " , round(errS,3),"; Percentage change on A: ",round(errA,3),".",sep=""))
    }
    A = A_new; S= S_new; theta = theta_new
    
    if(errA < espli1 & errS < espli2)
    {
      if(!silent){cat("Converaged!")}
      if(preprocess){ 
        A = Yraw %*% t(S) %*% solve(S%*%t(S)) 
      }else{
        A = Y %*% t(S) %*% solve(S%*%t(S)) 
      }
      return(list(Conver = TRUE, A=A,S=S,theta=theta))
    }
    Iter = Iter + 1
  }
  
  if(!silent){cat("Failed to converge!")}
  if(preprocess){ 
    A = Yraw %*% t(S) %*% solve(S%*%t(S)) 
  }else{
    A = Y %*% t(S) %*% solve(S%*%t(S)) 
  }
  return(list(Conver=FALSE, A=A, S=S, theta=theta))
}
