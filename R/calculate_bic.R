# Function to calculate BIC
calculate_bic = function(Y, A, S)
  # Y: connectivity data of dimension N x K, N is number of replicates/samples, K is number of edges. 
  # A: mixing matrix of dimension N x q, where q is the number of latent sources.
  # S: latent sources of dimension q by K. 
{
  N = dim(Y)[1]
  K = dim(Y)[2]
  
  # define a norm function
  norm_vec = function(x) sum(x^2)
  
  # calculate the residual standard deviation
  sigma = sqrt(1/(N*K)*sum(apply(Y - A%*%S, MARGIN = 1, norm_vec)))
  
  # calculate loglikelihood
  loglike = 0
  for (j in 1:N){
    mean = as.vector(A[j,] %*% S)
    loglike = loglike - 2 * sum(log(dnorm(Y[j,], mean, sigma)))
  }
  
  # calculate logN*sum(||S||0)
  L0 = log(N)*sum(abs(S)>1e-1)
  
  # calculate bic
  BIC = loglike + L0
  return(BIC)
}
