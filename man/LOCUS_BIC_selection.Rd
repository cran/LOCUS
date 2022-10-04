\name{LOCUS_BIC_selection}
\alias{LOCUS_BIC_selection}
\title{
BIC-based Hyper-parameters selection for LOCUS}
\description{
This function is to conduct the BIC-based hyper-parameters selection for LOCUS.
}
\usage{
LOCUS_BIC_selection(Y, q, V, MaxIteration=50, penalty="SCAD", 
phi_grid_search=seq(0.2, 1, 0.2), rho_grid_search=c(0.95), 
espli1=0.001, espli2=0.001, save_LOCUS_output=TRUE, 
preprocess=TRUE)
}

\arguments{
  \item{Y}{Group-level connectivity data from N subjects, which is of dimension N x p, where p is number of edges. Each row of Y represents a subject's vectorized connectivity matrix by \code{Ltrans} function. }
  \item{q}{Number of ICs/subnetworks to extract.}
  \item{V}{Number of nodes in the network. Note: p should be equal to V(V-1)/2. }
  \item{MaxIteration}{Maximum number of iteractions.}
  \item{penalty}{The penalization approach for uniform sparsity, which can be \code{NULL}, \code{SCAD}, \code{L1}, \code{Hardthreshold}. }
  \item{phi_grid_search}{Grid search candidates for tuning parameter of uniform sparse penalty.}
  \item{espli1}{Toleration for convergence on mixing coefficient matrix, i.e. A.}
  \item{espli2}{Toleration for convergence on latent sources, i.e. S.}
  \item{rho_grid_search}{Grid search candidates for tuning parameter for selecting number of ranks in each subnetwork's decomposition. }
  \item{save_LOCUS_output}{Whether to save LOCUS output from each grid search.}
  \item{preprocess}{Whether to preprocess the data, which reduces the data dimension to \code{q} and whiten the data.}
}

\details{
  In Wang, Y. and Guo, Y. (2023), the tuning parameters for learning the LOCUS model include \eqn{\phi, \rho}. The BIC-type criterion is proposed to select those parameters. 
  
\deqn{BIC =  -2 \sum_{i=1}^N log \{g(y_i; \sum_{l=1}^{q} \hat{a}_{il} \hat{s}_l, \hat{\sigma}^2 I_p)\} + log(N) \sum_{l=1}^{q}\|\hat{s}_l\|_0}

where \eqn{g} denotes the pdf of a multivariate Gaussian distribution, \eqn{\hat{\sigma}^2 = \frac{1}{Np}\sum_i \|y_i-\sum_{l=1}^{q}\hat{a}_{il}\hat{s}_{l}\|_2^2}, \eqn{\|\cdot\|_0} denotes the \eqn{L_0} norm . This criterion balances between model fitting and model sparsity. 
}
\value{
   \item{bic_tab}{BIC values per phi and rho.}
   \item{LOCUS_results}{LOCUS output, if save_LOCUS_output is TRUE. }
}
\examples{
## Simulated the data to use. 
V = 50
S1 = S2 = S3 = matrix(0,ncol = V,nrow = V)
S1[5:20,5:20] = 4;S1[23:37,23:37] = 3;S1[40:48,40:48] = 3
S2[15:20,] = -3;S2[,15:20] = -3
S3[15:25,36:45] = 3; S3[36:45,15:25] = 3
Struth = rbind(Ltrans(S1,FALSE) , Ltrans(S2,FALSE), Ltrans(S3,FALSE))
set.seed(100)
Atruth = matrix(rnorm(100*3),nrow=100,ncol=3)
Residual = matrix(rnorm(100*dim(Struth)[2]),nrow=100)
Yraw = Atruth\%*\%Struth + Residual

##### Run Locus on the data ##### 
Locus_bic_result = LOCUS_BIC_selection(Yraw,3,V)
print(Locus_bic_result$bic_tab)
# line plot
plot(Locus_bic_result$bic_tab[,2], Locus_bic_result$bic_tab[,3], type = "b",
     xlab = "phi", ylab = "BIC")
     
# visualize the best result based on BIC
idx = which.min(Locus_bic_result$bic_tab[,3])
oldpar = par(mfrow=c(2,3))
for(i in 1:3){image(Ltrinv(Struth[i,], V, FALSE))}
for(i in 1:3){image(Ltrinv(Locus_bic_result$LOCUS_results[[idx]]$LOCUS$S[i,], 
V, FALSE))}
par(oldpar)
}
