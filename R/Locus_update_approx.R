Locus_update_approx <- function(Y,A,theta,penalt = NULL,lambda_ch = 0.5,gamma = 3,imput_method = "Previous",silent = FALSE)
{
  # An extremely efficient approximation method with potentially higher performance. 
  if(is.null(penalt))
  {
    if(!silent)
      cat("Locus without penalty.")
  }else{
    if(!silent)
      cat(paste("Locus with", penalt,"penalty."))
  }
  
  theta_new = list()
  
  K = dim(Y)[2]
  V = (sqrt(1+8*K)+1)/2
  N = dim(Y)[1]
  q = length(theta)
  R = vector()
  for(curr_ic in 1:q)
  {
    theta_ic = theta[[curr_ic]]
    R[curr_ic] = dim(theta_ic$X_l)[1]
    S_lold = t(theta_ic$M_l%*%Y)
    
    if(is.null(penalt))
    {
      S_new_0 = S_lold
    }else if(penalt == "SCAD")
    {
      if(gamma<=2){warning("Gamma needs to be > 2!");gamma = 2.01}
      
      S_new_0= SCAD_func(S_lold,lambda_ch = lambda_ch  ,gamma = gamma)
      S_new_0 = S_new_0 /sd(S_new_0)*sd(S_lold)
      
    }else if(penalt == "Hardthreshold")
    {
      S_new_0 = S_lold*(abs(S_lold)>=lambda_ch)
    }else if(penalt == "L1")
    {
      S_new_0 = sign(S_lold)*(abs(S_lold)-lambda_ch)*(abs(S_lold)>=lambda_ch)
      S_new_0 = S_new_0 /sd(S_new_0)*sd(S_lold)
    }else
    {
      stop("No Penalty available!")
    }
    
    if(imput_method == "Previous"){
      Sl = Ltrinv(S_new_0,V,F) + diag(diag(t( theta_ic$X_l)%*%diag(theta_ic$lam_l)%*%theta_ic$X_l ))
    }else if(imput_method == "Average"){ 
      Sl = Ltrinv(S_new_0,V,F) + diag( rep(mean(S_new_0),V )) 
    }else{
      stop("No Imputation available!")
    }
    eigenSl = eigen(Sl)
    orderEigen = order(abs(eigenSl$values),decreasing = T)
    Rl = R[curr_ic]
    eigenset = orderEigen[1:Rl]
    
    theta_new[[curr_ic]] = list()
    theta_new[[curr_ic]]$lam_l = eigenSl$values[ eigenset ]
    if( theta_new[[curr_ic]]$lam_l[1]<0 ) 
    {
      theta_new[[curr_ic]]$lam_l = -1*theta_new[[curr_ic]]$lam_l
    }
    for(j in 1:Rl)
    {
      theta_ic$X_l[j,]= eigenSl$vectors[,eigenset[j]]
    }
    theta_new[[curr_ic]]$X_l = theta_ic$X_l
  }
  
  # Update A 
  l = 1; Snew = array(dim=c(q,K))
  while(l <= q)
  {
    Snew[l,] = Ltrans(t(theta_new[[l]]$X_l)%*% diag(theta_new[[l]]$lam_l) %*%theta_new[[l]]$X_l,F)  # K x M
    l = l+1
  }
  
  ## estimate A
  Anew = Y%*%t(Snew) %*% solve(Snew%*%t(Snew)) 
  for(i in 1:q)
  {
    ai = sd(Anew[,i])
    theta_new[[i]]$lam_l = theta_new[[i]]$lam_l * ai
    Anew[,i] =Anew[,i] / ai
    Snew[i,] = Snew[i,] * ai
  }
  
  # Save m_l, X_l into theta2_new:
  Mnew =  solve(t(Anew)%*%Anew)%*%t(Anew)
  for(l in 1:q)
  {
    theta_new[[l]]$M_l = Mnew[l,]
  }
  return(list(A = Anew,S=Snew,theta = theta_new))
}
