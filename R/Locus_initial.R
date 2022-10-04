Locus_initial <- function(Y,q,V,rho=0.95,R = NULL,maxIter = 100)
{ 
  ICcorr = icaimax(t(Y),nc = q,center=FALSE,maxit=maxIter)
  S_ini = matrix(0,ncol=dim(ICcorr$S)[1],nrow=q)
  theta_ini = list()
  for(i in 1:q)
  {
    Sl = Ltrinv( ICcorr$S[,i],V,FALSE)
    Sl = Sl + diag( rep(mean(ICcorr$S[,i]),V ))
    eigenSl = eigen(Sl)
    orderEigen = order(abs(eigenSl$values),decreasing = TRUE)
    if(is.null(R))
    {
      Rl = 2
      while( TRUE )
      {
        eigenset = orderEigen[1:Rl]
        imgeRL = eigenSl$vectors[,eigenset]%*% diag(eigenSl$values[eigenset])%*% t(eigenSl$vectors[,eigenset])
        # image( imgeRL ) 
        if(cor(Ltrans(imgeRL,FALSE),ICcorr$S[,i]) > rho) break
        Rl = Rl + 1
      }
    }else
    {
      Rl = R[i]; eigenset = orderEigen[1:Rl]
    }
    theta_ini[[i]] = list()
    
    theta_ini[[i]]$lam_l = eigenSl$values[ eigenset ]
    if( theta_ini[[i]]$lam_l[1]<0 ){theta_ini[[i]]$lam_l = -1*theta_ini[[i]]$lam_l}    
    theta_ini[[i]]$X_l = matrix(0,ncol = V, nrow = Rl)
    for(j in 1:Rl)
    {
      theta_ini[[i]]$X_l[j,] = eigenSl$vectors[,eigenset[j]]
    }
    S_ini[i,] = Ltrans( t(theta_ini[[i]]$X_l)%*%diag(theta_ini[[i]]$lam_l)%*%theta_ini[[i]]$X_l,FALSE)
  }
  
  # compute initial value of A
  A_ini = Y%*%t(S_ini)%*%solve(S_ini%*%t(S_ini))
  
  # scale up
  for(l in 1:q){
    # unit norm each column of A
    # scaleL = sqrt(sum(A_ini[,l]^2))
    scaleL = sd(A_ini[,l])
    A_ini[,l] = A_ini[,l] / scaleL
    S_ini[l,] = S_ini[l,] * scaleL
    # scale X_l correspondingly
    theta_ini[[l]]$X_l = theta_ini[[l]]$X_l * sqrt(scaleL)
  }
  
  # since after preprocessing, A_tilde is orthogonal
  # compute A_tilde transpose/inverse (g-inverse of A-tilde)
  # if X has full column rank, (X'X)^(-1)X' is its g-inverse
  # why use g-inverse here: since A_ini is N*q
  # afterwards, in the update process, we actually can just use A_tilde transpose(after scale), or g-inverse
  M_ini =  solve(t(A_ini)%*%A_ini)%*%t(A_ini) # g-inverse of A
  for(l in 1:q){
    theta_ini[[l]]$M_l = M_ini[l,]
  }
  
  return(list(A=A_ini,theta = theta_ini,S = S_ini))
}

