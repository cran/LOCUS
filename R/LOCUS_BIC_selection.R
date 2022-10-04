# Use BIC method to select the hyper-parameters of phi and rho. 
LOCUS_BIC_selection <- function(Y, q, V, MaxIteration=50, penalty="SCAD", phi_grid_search=seq(0.2, 1, 0.2), 
                      rho_grid_search=c(0.95), espli1=0.001, espli2=0.001, save_LOCUS_output=TRUE, 
                      preprocess=TRUE)
  # Y: connectivity data of dimension N x K, N is number of subjects, K is number of edges. 
  # q: Number of subnetworks to extract.
  # V: Number of nodes in network.
  # MaxIteration = 100: number of maximum iterations.
  # espli1=0.001, espli2=0.001: toleration for convergence on S and A. 
  # penalty = "SCAD": sparsity penalty, this can be NULL, SCAD, or L1.
  # phi_grid_search: grid search candidates of tuning parameter for penalty
  # rho_grid_search: grid search candidates of tuning parameter for selecting number of ranks in each subnetwork's decomposition. 
{
  
  # demean or preprocess Y to speed up 
  Y = sweep(Y,2,apply(Y,2,mean),"-") 
  if(preprocess){
    Y = Locus_preprocess(Y, q)
  }
  
  # run grid search  
  LOCUS_results = list()
  bic_tab = matrix(0, ncol = 3, nrow=length(rho_grid_search) * length(phi_grid_search))
  k = 1
  for(rho in rho_grid_search)
  {
    for(phi in phi_grid_search)
    {
      message('Running LOCUS for rho=', rho, ', phi=', phi, '...')
      LOCUS_output = LOCUS_internal_in_bic(Y, q=q, V=V, MaxIteration=MaxIteration, penalty=penalty, phi=phi, approximation=TRUE,
                           preprocess=FALSE, espli1=espli1, espli2=espli1, rho=rho, silent=TRUE, demean=FALSE)
      bic_value = calculate_bic(Y, LOCUS_output$A, LOCUS_output$S)
      bic_tab[k,] = c(rho, phi, bic_value)
      if(save_LOCUS_output)
      {
        LOCUS_results[[k]] = list(LOCUS=LOCUS_output, rho=rho, phi=phi)
      }
      k = k+1
      
    }
  }
  
  if(save_LOCUS_output)
  {
    colnames(bic_tab) = c("rho", "phi", "bic_value")
    return(list(bic_tab = bic_tab, LOCUS_results = LOCUS_results))
  }
  return(list(bic_tab = bic_tab))
}
