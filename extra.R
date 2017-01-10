
##########################################################
# Step 6: Run calibration with diagnostics
##########################################################
computeModelZCBPricesWithDiagnostics <- function(eigenvals)
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
  
  B.tT <- B.tau %*% x.0
  print("B(t,T)^T*x_t:")
  print(B.tT)
  print('######################')
  
  A.tT <- c(constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[1]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[2]),
            constructAffineAMatrix(eigenvals, kappa.Q, F.matrix, maturities[3]))
  print("A(t,T):")
  print(A.tT)
  print('######################')
  
  model.zero.coupon.bond.prices <- exp(A.tT + as.vector(B.tT))
  
  print("Model Zero Coupon Bond Prices:")
  print(model.zero.coupon.bond.prices)
  print('######################')
  
  return (model.zero.coupon.bond.prices)
}

calibrateKappaQWithDiagnostics <- function()
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
  
  model.zcb.prices <- computeModelZCBPricesWithDiagnostics(kappa.eigenvals)
  
  model.yields <- fromPricesToYields(model.zcb.prices, maturities)
  
  print("Model Yields:")
  print(model.yields)
  print('######################')
  
}




