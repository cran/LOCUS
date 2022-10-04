Locus_update <- function(Y,A,theta,penalt = NULL,lambda_ch = 0.5,gamma = 3,silent = FALSE)
{
  # Update latent channels X's with flat prior 
  # A denotes mixing matrix
  # theta contains IC-specific parameters and channels
  if(is.null(penalt))
  {
    if(!silent)
      cat("Locus without penalty.")
  }else{
    if(!silent)
      cat(paste("Locus with", penalt,"penalty."))
  }
  
  theta_new = list()
  #require(MASS)
  
  K = dim(Y)[2]         # Number of edges
  V = (sqrt(1+8*K)+1)/2
  N = dim(Y)[1]
  q = length(theta)
  
  rmat = matrix(rep(1:V,V),ncol=V)
  cmat = t(rmat)
  Lcoorsym = matrix(0,ncol=2,nrow=V*(V-1)/2)
  Lcoorsym[,1] = Ltrans(rmat,F)
  Lcoorsym[,2] = Ltrans(cmat,F)
  
  ## For each IC, conditioned on others, estimate latent channels X: 
  
  newS = array(dim=c(q,K))
  for(curr_ic in 1:q)
  {
    theta_new[[curr_ic]] = list()
    
    # Update X: 
    theta_ic = theta[[curr_ic]]
    Dinverse = 1/theta[[curr_ic]]$lam_l           # D^(-1): of length R
    theta[[curr_ic]]$X_l[is.na(theta[[curr_ic]]$X_l)] = 0
    Yic = t(theta[[curr_ic]]$M_l%*%Y)             # p x 1
    v = 1 
    while(v <= V)
    {
      Hlv = t(theta[[curr_ic]]$X_l[,-v])                    # X(-v) which is V-1 x R
      yvpen = Yic[which(Lcoorsym[,1] == v | Lcoorsym[, 2] == v),]               # Y 
      Sigmalv = psdmat_inverse( t(Hlv)%*%Hlv )
      if(is.null(penalt)) 
      {
        beta = yvpen
      }else if(penalt == "SCAD")
      {
        if(gamma<=2){warning("Gamma needs to be > 2!");gamma = 2.01}
        beta= SCAD_func(yvpen,lambda_ch = lambda_ch  ,gamma = gamma)
      }else if(penalt == "L1")
      {
        beta = sign(yvpen)*(abs(yvpen)-lambda_ch)*(abs(yvpen)>=lambda_ch)
      }else if(penalt == "Hardthreshold")
      {
        beta = yvpen*(abs(yvpen)>=lambda_ch)
      }else
      {
        stop("No Penalty available!")
      }
      if(sd(beta) == 0)
      {
        theta[[curr_ic]]$X_l[,v] = 0
      }else
      {
        theta[[curr_ic]]$X_l[,v] = (Dinverse*(Sigmalv%*%t(Hlv)%*%beta)) /sd(beta)*sd(yvpen) # Rx1
      }
      v = v+1
    }
    theta_new[[curr_ic]]$X_l = theta[[curr_ic]]$X_l
    
    # Update D: 
    Xstarstack = t(apply(theta[[curr_ic]]$X_l,1,function(x){ x = matrix(x,ncol=1);return(Ltrans(x%*%t(x),F)) }))
    
    if(is.null(penalt)) 
    {
      beta = Yic
    }else if(penalt == "SCAD")
    {
      if(gamma<=2){warning("Gamma needs to be > 2!");gamma = 2.01}
      beta= SCAD_func(Yic,lambda_ch = lambda_ch  ,gamma = gamma)
    }else if(penalt == "L1")
    {
      beta = sign(Yic)*(abs(Yic)-lambda_ch)*(abs(Yic)>=lambda_ch)
    }else if(penalt == "Hardthreshold")
    {
      beta = Yic*(abs(Yic)>=lambda_ch)
    }
    
    theta_new[[curr_ic]]$lam_l = as.numeric(solve(Xstarstack%*%t(Xstarstack))%*%Xstarstack%*%beta /sd(beta)*sd(Yic)) # Rx1
    
    newS[curr_ic,] = Ltrans(t(theta_new[[curr_ic]]$X_l)%*% 
                              diag(theta_new[[curr_ic]]$lam_l) %*%
                              theta_new[[curr_ic]]$X_l,F)  
  }
  
  # Update A
  newA = Y%*%t(newS) %*% solve(newS%*%t(newS)) 

  for(i in 1:q)
  {
    ai = sd(newA[,i])
    theta_new[[i]]$lam_l = theta_new[[i]]$lam_l * ai
    newA[,i] = newA[,i] / ai
    newS[i,] = newS[i,] * ai
  }
  
  # Save m_l, X_l into theta_new:
  Mnew =  solve(t(newA)%*%newA)%*%t(newA) # Generalized inverse of A
  for(l in 1:q)
  {
    theta_new[[l]]$M_l = Mnew[l,]
  }
  return(list(A=newA,theta=theta_new,S=newS))
}
