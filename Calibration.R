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
u.r <- 0.05

# (3) Mean Reversion Levels
mean.rev.levels <- matrix(c(0.01, 0.02, 0.03), nrow=3, ncol=1)

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
constructKappaQMatrix <- function(eigenvals, mats)
{
  F.matrix <- constructFMatrix(eigenvals, mats)
  kappa.Q <- beta.tilde.inv %*% F.matrix %*% eigenvals %*% solve(F.matrix) %*% beta.tilde
  return (kappa.Q)
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
  L.inv <- solve(L)
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
  I3b <- -0.5 * t(m.vector) %*% kappa.Q.inv %*% C %*% t(solve(L)) %*% D.T %*% solve(eigenvals) %*% t(L) %*% m.vector
  I3c <- -0.5 * t(m.vector) %*% kappa.Q.inv %*% L %*% D.T %*% L.inv %*% C %*% solve(t(kappa.Q)) %*% m.vector
  I3d <- 0.5 * tau * t(m.vector) %*% kappa.Q.inv %*% C %*% solve(t(kappa.Q)) %*% m.vector
  return (I1 + I2 + I3a + I3b + I3c + I3d)
}

computeModelZCBPrices <- function(eigenvals)
{
  kappa.Q <- constructKappaQMatrix(eigenvals, maturities)
  print("Kappa matrix in Q measure:")
  print(kappa.Q)
  print('######################')
  
  m.vector <- computeMVector()
  print("m vector")
  print(m.vector)
  print('######################')
  
  F.matrix <- constructFMatrix(eigenvals, maturities)
  print("F Matrix:")
  print(F.matrix)
  print('######################')
  
  B.tau <- rbind(constructAffineBMatrix(kappa.Q, maturities[1]),
                 constructAffineBMatrix(kappa.Q, maturities[2]),
                 constructAffineBMatrix(kappa.Q, maturities[3]))
  print("B(t,T)^T Matrix:")
  print(B.tau)
  print('######################')
  
  #Get the state variables for t=0
  x.0 <- c(tail(t(two.pca$PrincipalComponents))[1,], 0.01)
  print("State variables for t=0:")
  print(x.0)
  print('######################')
  
  B.x.0 <- B.tau %*% x.0
  print("B(t,T)^T*x_t:")
  print(B.x.0)
  print('######################')
  
  A.tT <- c(constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[1]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[2]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[3]))
  print("A(t,T):")
  print(A.tT)
  print('######################')
  
  model.zero.coupon.bond.prices <- exp(A.tT + as.vector(B.x.0))
  
  print("Model Zero Coupon Bond Prices:")
  print(model.zero.coupon.bond.prices)
  print('######################')
  
  return (model.zero.coupon.bond.prices)
}

calibrateKappaQ <- function()
{
  print('######################')
  #Get the initial eigenvalues of kappa
  kappa.eigenvals <- kappa.eigenvals.initial
  print("Initial Eigenvalues:")
  print(kappa.eigenvals)
  print('######################')
  
  #Get the initial term structure. We assume that the valuation date or calibration date
  #is the last date on which the data is available.
  initial.term.structure <- tail(three.yields, 1)
  print('Initial Term Structure:')
  print(initial.term.structure)
  print('######################')
  
  initial.zero.coupon.bond.prices <- exp(-initial.term.structure*maturities)
  print('Initial ZCB Prices:')
  print(initial.zero.coupon.bond.prices)
  print('######################')
  
  model.zcb.prices <- computeModelZCBPrices(kappa.eigenvals)
  
  model.yields <- fromPricesToYields(model.zcb.prices, maturities)
  
  print("Model Yields:")
  print(model.yields)
  print('######################')
  
}

calibrateKappaQ()


