Locus_preprocess <-function(Ynew,q) 
{
  N = dim(Ynew)[1]
  eigenA = eigen( Ynew %*% t(Ynew), T )
  whitenmat = diag( (eigenA$values[1:q] - mean(eigenA$values[(q+1):N]))^(-0.5)  ) %*% t(eigenA$vectors[,1:q])
  Ynew = whitenmat%*%Ynew # the whitened A matrix
  Ynew = Ynew/sd(Ynew)*5
  return(Ynew)
}
