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

##########################################################
# Step 1: Do the PCA on N (=2) yields and N+1 (=3) yields
# Step 2: Compute the Beta matrix
# Both these steps are done in the file BetaComputation.r
##########################################################
source(paste(working_dir,"BetaComputation.R",sep=""))

##########################################################
# Step 3: Calibrate for the kappa matrix under the Q measure
##########################################################
#Initialize eigenvalues of Kappa
kappa.eigenvals.initial <- diag(c(0.037,0.0631,0.0630))

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
  kappa.Q <- beta.tilde %*% F.matrix %*% eigenvals %*% solve(F.matrix) %*% beta.tilde.inv
  return (kappa.Q)
}

#Given the 3 eigenvalues (l1,l2,l3) and their maturities (tau1, tau2, tau3)
#We build the 3 x 1 g vector
constructGVector <- function(eigenvals, mats)
{
  #Get the matrix of eigenvectors for 3 yields
  omega_tilde <- three.pca$Omega
  #Construct the F matrix using the eigenvalues and maturities provided
  F.matrix <- constructFMatrix(eigenvals, mats)
  #Unity 3 x 1 matrix
  unity = matrix(1, nrow=3, ncol=1)
  #Return the g vector
  return (t(omega_tilde) %*% solve(t(F.matrix)) %*% unity)
}

#Given the 3 eigenvalues (l1,l2,l3) and their maturities (tau1, tau2, tau3)
#We build the 3 x 1 g vector
constructBtTMatrixForAffineFunction <- function(kappa.Q, g.vector, tau)
{
  return ((expm(-t(kappa.Q)*tau) - diag(x=1,nrow=3,ncol=3)) %*% solve(t(kappa.Q)) %*% g.vector)
}

constructAtTMatrixForAffineFunction <- function(kappa.Q, g.vector, tau)
{
  return (c(0,0,0))
}

computeModelZCBPrices <- function(eigenvals)
{
  kappa.Q <- constructKappaQMatrix(eigenvals, maturities)
  print("Kappa matrix in Q measure:")
  print(kappa.Q)
  print('######################')
  
  g.vector <- constructGVector(eigenvals, maturities)
  print("g vector")
  print(g.vector)
  print('######################')
  
  F.matrix <- constructFMatrix(eigenvals, maturities)
  print("F Matrix:")
  print(F.matrix)
  print('######################')
  
  B.tau <- cbind(constructBtTMatrixForAffineFunction(kappa.Q, g.vector, maturities[1]),
                 constructBtTMatrixForAffineFunction(kappa.Q, g.vector, maturities[2]),
                 constructBtTMatrixForAffineFunction(kappa.Q, g.vector, maturities[3]))
  print("B(t,T) Matrix:")
  print(B.tau)
  print('######################')
  
  #Get the state variables for t=0
  x.0 <- tail(t(three.pca$PrincipalComponents))[1,]
  print("State variables for t=0:")
  print(x.0)
  print('######################')
  
  B.x.0 <- t(B.tau) %*% x.0
  print("B(t,T)*x_t:")
  print(B.x.0)
  print('######################')
  
  A.tau <- constructAtTMatrixForAffineFunction(kappa.Q, g.vector, maturities)
  model.zero.coupon.bond.prices <- exp(A.tau + as.vector(B.x.0))
  
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
  initial.term.structure <- tail(three.yields)[1]
  print('Initial Term Structure:')
  print(initial.term.structure)
  print('######################')
  
  initial.zero.coupon.bond.prices <- exp(-initial.term.structure*maturities)
  print('Initial ZCB Prices:')
  print(initial.zero.coupon.bond.prices)
  print('######################')
  
  model.zcb.prices <- computeModelZCBPrices(kappa.eigenvals)
}

calibrateKappaQ()


