
load_yields_data <- function(data_file)
{
  #Read the yield curve data in the CSV format, remove the NA values
  yields <- na.omit(read.csv(data_file))
  #Convert the yields into xts
  #yields.ts <- xts(yields[,-1], order.by=as.Date(as.character.Date(yields[,1]),format="%d/%m/%Y"))
  yields.ts <- xts(yields[,-1], order.by=as.Date(yields[,1]))
  
  #Get the list of maturities from the yield data
  maturities <- as.numeric(substring(colnames(yields.ts),6,7))

  if(use.plotly)
  {
    #Plot the yield surface
    yield.curves.plot <- plot_ly(z=as.matrix(yields[,-1])) %>% add_surface
    yield.curves.link <- plotly_POST(yield.curves.plot, filename="RawYieldCurve")
  }
  
  return(yields.ts)
}

