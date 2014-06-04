PC.Analysis = function(table)  
{
  #read the table
  x <- table
  #First we have to scale the table by: 
  #x-scaled(ij) = (x - column mean)/column range 
  #the centering is automatically done, but we have to find the range.
  #First we have to find the number of columns, because that determines
  #the number of values that we have to put into the scale.vector:
  xncol <- ncol(x)
  #Now we find the scaling values, which are basically the column ranges:
  #First we create an empty vector:
  #range.vector=c()
  #Find the range values and place in a vector:
  for (j in 1:xncol){
    min.col [j] <- min(x [,j])
    max.col [j] <- max(x [,j])    
  }
  rangediff <- max.col - min.col
  y <- scale(x, center = TRUE, scale=rangediff)
  ysvd <- svd(y)
  
  #Calculate squared of the singular value...
  sqr.ysvd <- ysvd$d^2
  #...and its %
  total.sqr.ysvd <- sum(sqr.ysvd)
  percent.ysvd <- (sqr.ysvd/total.sqr.ysvd)*100
  
  #Plot the vector coordinates
  plot <- plot(ysvd$v)
  
  ds = sum(y * y) #Data scatter
  mu1 <- round(ysvd$d[1] ^2) #Mu
  contr <- (mu1 / ds)*100 #Contribution of first component
  
  answer <- list(ysvd$d, sqr.ysvd, percent.ysvd, ysvd$u, ysvd$v, ds, mu1, contr, plot)
  names(answer)[[1]] <- "Singular.value"
  names(answer)[[2]] <- "Singular.value.squared"
  names(answer)[[3]] <- "Singular.value.percent"
  names(answer)[[4]] <- "Scatter.coordinates"
  names(answer)[[5]] <- "PC.vector.coordinates"
  names(answer)[[6]] <- "Data.scatter"
  names(answer)[[7]] <- "Mu"
  names(answer)[[8]] <- "First.component.Contribution.%"
  names(answer)[[9]] <- "Plot"  
  return(answer)
}