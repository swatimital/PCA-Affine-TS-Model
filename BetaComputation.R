#Do the PCA Analysis on time series data
pca_analysis <- function(Raw_Data)
{
  mu <- as.vector(apply(Raw_Data, 2, mean))
  X <- t(t(Raw_Data)-mu)
  Cov_Matrix <- cov(X)
  aux <- eigen(Cov_Matrix)
  Omega <- aux$vectors
  PrincipalComponents <- t(Omega) %*% t(X)
  Reconstr_Data <- t(Omega %*% PrincipalComponents + mu)
  return (list(Reconstr_Data = Reconstr_Data, Omega = Omega, PrincipalComponents = PrincipalComponents, D = diag(aux$values)))
}

#Do the PCA on all yields and plot the first 3 principal components
pca.on.yields <- pca_analysis(yields.ts)
matplot(x=1:length(pca.on.yields$Omega[,1]),y=pca.on.yields$Omega[,1:3],type='l',col=c('black','red','blue'))

#This will compute the state variables using N=2 yields, we have chosen 2y, 5y
two.pca <- pca_analysis(two.yields)
x <- two.pca$PrincipalComponents
dx <- t(diff(t(x)))

#This will compute the state variables using N+1=3 yeilds, we have chosen 2y, 5y and 10y
three.pca <- pca_analysis(three.yields)
x.tilde <- three.pca$PrincipalComponents
dx.tilde <- t(diff(t(x.tilde)))

#Get the covariance of the 3 yields
sigma_p <- cov(three.yields)
#Get the matrix of eigenvectors for 3 yields
omega_tilde <- three.pca$Omega
#Compute dx.dp^T
dx.dp <- dx.tilde %*% t(dx)
dx.dp <- cbind(dx.dp,c(0,0.5,0))
#Compute beta and beta inverse matrices
beta.tilde <- omega_tilde*dx.dp*solve(sigma_p)
beta.tilde.inv <- solve(beta.tilde)