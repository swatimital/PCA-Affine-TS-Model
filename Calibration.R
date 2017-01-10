working_dir <- "/Users/swatimital/GitHub/Thesis/TermStructure/AffineModel/"
source(paste(working_dir,"Settings.r",sep=""))
source(paste(working_dir,"DataLoader.r",sep=""))

#Fed Data filename
fed.data.filename <- paste(working_dir, 'FED-SVENY.csv', sep="")
#Select 2y,5y,10y yields
maturities <- c(2,5,10)
#Load the Zero Coupon Yields
yields.ts <- load_yields_data(fed.data.filename)/100
#These are the 3 yield curves that are selected (corresponding to y tilde in the paper)
three.yields <- cbind(yields.ts$SVENY02, yields.ts$SVENY05, yields.ts$SVENY10)
#Compute the corresponding zero coupon bond prices
three.zcb.prices <- exp(-three.yields*maturities)
#These are the 2 yield curves that are selected (corresponding to y in the paper)
two.yields <- cbind(yields.ts$SVENY02, yields.ts$SVENY05)
#Compute the corresponding zero coupon bond prices
two.zcb.prices <- exp(-two.yields*maturities[1:2])
#Short Rate at t=T, t=T-1, t=T-2
r.t <- as.matrix(c(0.0043678, 0.0043733, 0.0043456), nrow=3, ncol=1)

##########################################################
# Step 1: Do the PCA on N (=2) yields and N+1 (=3) yields
# Step 2: Compute the Beta matrix
# Both these steps are done in the file BetaComputation.r
##########################################################
source(paste(working_dir,"BetaComputation.R",sep=""))

##########################################################
#Step 3: Calibration User Inputs
##########################################################
# (1) Initialize eigenvalues of Kappa
kappa.eigenvals.initial <- diag(c(0.01,0.5,0.6))
#kappa.eigenvals.initial <- diag(c(0.037,0.631,0.630))

# (2) u.r
u.r <- 0.01

# (3) Mean Reversion Levels
mean.rev.levels <- matrix(c(0.005, 0.008, 0.01), nrow=3, ncol=1)

# (4) 3 most recent values of the MPR in order to calibrate the
# m vector
lambda.t <- c(-0.001, -0.0012, -0.0013)
##########################################################


##########################################################
# Step 4: Define Helper Functions that are not dependent 
# on eigenvalues
##########################################################
computeEtaVector <- function()
{
  return (matrix(as.vector(apply(three.yields, 2, mean)), nrow=3, ncol=1))
}

computeMVector <- function()
{
  x.t <- tail(t(two.pca$PrincipalComponents),3)
  x.t <- x.t[nrow(x.t):1,]
  rhs.vec <- r.t - rep(u.r,3)
  x.t <- cbind(x.t, lambda.t)
  return (solve(x.t) %*% rhs.vec)
}

computeSpMatrix <- function()
{
  return (beta.tilde.inv %*% three.pca$Omega %*% sqrt(three.pca$EigenValues))
}

computeThetaVector <- function()
{
  return (beta.tilde.inv %*% (three.pca$Omega %*% mean.rev.levels + computeEtaVector()))
}

fromYieldToPrices <- function(yields, mats)
{
  return (exp(-yields*mats))
}

fromPricesToYields <- function(prices, mats)
{
  return (-log(prices)/mats)
}

##########################################################
# Step 5: Calibrate for the kappa matrix under the Q measure
##########################################################

#Given the 3 eigenvalues (l1,l2,l3) and their maturities (tau1, tau2, tau3)
#We build the 3 x 3 F matrix
constructFMatrix <- function(eigenvals, mats)
{
  F <- matrix(0,nrow=3,ncol=3)
  for (i in 1:3)
  {
    for(j in 1:3)
    {
      F[i,j] <- (1-exp(-eigenvals[j,j]*mats[i]))/(eigenvals[j,j]*mats[i])
    }
  }
  return(F)
}

#Given the 3 eigenvalues (l1,l2,l3) and their maturities (tau1, tau2, tau3)
#We build the 3 x 3 mean reversion matrix K under Q measure 
constructKappaQMatrix <- function(eigenvals, f.matrix, mats)
{
  return (beta.tilde.inv %*% f.matrix %*% eigenvals %*% solve(f.matrix, tol=1e-19) %*% beta.tilde)
}

#Given the 3 eigenvalues (l1,l2,l3) and their maturities (tau1, tau2, tau3)
#We build the 3 x 1 g vector
constructAffineBMatrix <- function(kappa.Q, tau)
{
  #return ((expm(-t(kappa.Q)*tau) - diag(x=1,nrow=3,ncol=3)) %*% solve(t(kappa.Q)) %*% m.vector)
  return (t(computeMVector()) %*% solve(kappa.Q) %*% (expm(-kappa.Q*tau)-diag(x=1,nrow=3,ncol=3)))
}

constructAffineAMatrix <- function(eigenvals, kappa.Q, F.matrix, tau)
{
  #Compute the intermediate values
  m.vector <- computeMVector()
  C <- computeSpMatrix()
  L <- beta.tilde.inv %*% F.matrix
  L.inv <- solve(L, tol=1e-19)
  D.T <- diag(x=as.vector(diag2vec(1-exp(-eigenvals*tau))/diag2vec(eigenvals)), nrow=3, ncol=3)
  M <- L.inv %*% C %*% t(L.inv)
  F.T <- matrix(0,nrow=3,ncol=3)
  kappa.Q.inv <- solve(kappa.Q)
  theta.p <- computeThetaVector()
  
  for (i in 1:3)
  {
    for (j in 1:3)
    {
      F.T[i,j] <- (M[i,j] * (1-exp(-(eigenvals[i,i]+eigenvals[j,j])*T))) / (eigenvals[i,i]+eigenvals[j,j])
    }
  }
  #A(t,T) = I1 + I2 + I3a + I3b + I3c + I3d
  I1 <- -u.r*tau
  I2 <- t(m.vector) %*% (L %*% solve(eigenvals) %*% D.T %*% eigenvals %*% L.inv %*% theta.p - theta.p * tau)
  I3a <- 0.5 * t(m.vector) %*% kappa.Q.inv %*% L %*% F.T %*% solve(eigenvals) %*% t(L) %*% m.vector
  I3b <- -0.5 * t(m.vector) %*% kappa.Q.inv %*% C %*% t(L.inv) %*% D.T %*% solve(eigenvals) %*% t(L) %*% m.vector
  I3c <- -0.5 * t(m.vector) %*% kappa.Q.inv %*% L %*% D.T %*% L.inv %*% C %*% solve(t(kappa.Q)) %*% m.vector
  I3d <- 0.5 * tau * t(m.vector) %*% kappa.Q.inv %*% C %*% solve(t(kappa.Q)) %*% m.vector
  return (I1 + I2 + I3a + I3b + I3c + I3d)
}

computeModelZCBPrices <- function(eigenvals, mats)
{
  F.matrix <- constructFMatrix(eigenvals, mats)
  kappa.Q <- constructKappaQMatrix(eigenvals, F.matrix, mats)
  B.tau <- constructAffineBMatrix(kappa.Q, mats[1])
  for (m in mats[-1])
  B.tau <- rbind(B.tau, constructAffineBMatrix(kappa.Q, m))
  #Get the state variables for t=0
  x.0 <- c(tail(t(two.pca$PrincipalComponents))[1,], lambda.t[1])
  B.tT <- B.tau %*% x.0
  A.tT <- c(constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[1]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[2]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[3]))
  model.zero.coupon.bond.prices <- exp(A.tT + as.vector(B.tT))
  return (model.zero.coupon.bond.prices)
}

min.RSS <- function(data, par)
{
  #Maturities of the market bond prices
  mats <- data$x
  #Convert the initial eigenvalues into diagonal matrix
  eigenvals <- vec2diag(par)
  F.matrix <- constructFMatrix(eigenvals, mats)
  kappa.Q <- constructKappaQMatrix(eigenvals, F.matrix, mats)
  
  B.tau <- rbind(constructAffineBMatrix(kappa.Q, mats[1]),
                 constructAffineBMatrix(kappa.Q, mats[2]),
                 constructAffineBMatrix(kappa.Q, mats[3]))
  
  #Get the state variables for t=0
  x.0 <- c(tail(t(two.pca$PrincipalComponents))[1,], lambda.t[1])
  B.tT <- B.tau %*% x.0
  A.tT <- c(constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, mats[1]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, mats[2]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, mats[3]))
  model.zero.coupon.bond.prices <- exp(A.tT + as.vector(B.tT))
  print( sum((model.zero.coupon.bond.prices - data$y)^2) )
  print(model.zero.coupon.bond.prices)
  return (sum((model.zero.coupon.bond.prices - data$y)^2))
}


#Get the initial eigenvalues of kappa
kappa.eigenvals <- kappa.eigenvals.initial
#Get the initial term structure. We assume that the valuation date or calibration date
#is the last date on which the data is available.
current.term.structure <- tail(three.yields, 1)
current.zcb.prices <- exp(-current.term.structure*maturities)

dat <- data.frame(x=maturities, y=t(current.zcb.prices))
colnames(dat) <- c("x","y")
fit <- optim(par=diag2vec(kappa.eigenvals.initial), min.RSS, data=dat, hessian=TRUE)
zcb.prices <- computeModelZCBPrices(vec2diag(fit$par), maturities)
zcb.prices

#Compute confidence intervals for l1,l2,l3
fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
interval<-data.frame(value=fit$par, upper=upper, lower=lower)
interval


