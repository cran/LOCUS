psdmat_inverse<-function(mat)
{
  # Originally defined for PSD matrices
  p = dim(mat)[1]
  eigendecomp = eigen(mat)
  
  #if(-min(eigendecomp$values)> 0.0001){print("Error: Matrix has negative eigenvalue!");stop()}
  
  if( min(eigendecomp$values)<=0.0001 )
  {
    message("Matrix has nearly zero eigenvalue, perturbation is added. ")
    # print(round(eigendecomp$values,3))
    perturb <- max(max(eigendecomp$values) - p * min(eigendecomp$values), 0)/(p - 1)
  }else
  {
    perturb = 0
  }

  return( (eigendecomp$vectors)%*%diag(1/(perturb+eigendecomp$values))%*%t(eigendecomp$vectors) )
}
