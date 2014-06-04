PC.Analysis.NN = function(table)  
{
  #read the table
  x <- table
  #First we have to scale the table by: 
  #x-scaled(ij) = (x - column mean)/column range 
  #the centering is automatically done, but we have to find the range.
  #First we have to find the number of columns, because that determines
  #the number of values that we have to put into the scale.vector:
  #xncol <- ncol(x)
  #Now we find the scaling values, which are basically the column ranges:
  #First we create an empty vector:
  #range.vector=c()
  #Find the range values and place in a vector:
  #for (j in 1:xncol){
  #  min.col [j] <- min(x [,j])
  #  max.col [j] <- max(x [,j])    
  #}
  #rangediff <- max.col - min.col
  y <- scale(x, scale=FALSE) 
  ysvd <- svd(x)
  
  #Calculate squared of the singular value...
  sqr.ysvd <- ysvd$d^2
  #...and its %
  total.sqr.ysvd <- sum(sqr.ysvd)
  percent.ysvd <- (sqr.ysvd/total.sqr.ysvd)*100
  #Calculate square root of the singular value...
  sqrt.ysvd <- sqrt(ysvd$d)
  
  #Find the contributions to the hidden factor Feature Loadings
  feature.loading <- ysvd$v
  flncol <- ncol(feature.loading)
  for (j in 1:flncol){
    feature.loading[,j] <- ysvd$v[,j] * -sqrt.ysvd[j]
  }
  #Add the column names
  PC.FeatureLoad <- list(NULL, paste("Feature.Loading.", 1:flncol, sep=""))
  dimnames(feature.loading) <- PC.FeatureLoad
  #Add the row names
  flnrow <- nrow(feature.loading)
  rownames(feature.loading) <- rownames(feature.loading, do.NULL=FALSE, prefix="C.")
  #For convenience, transpose the matrix
  transp.feature.loading <- t(feature.loading)
  
  #Find the Factor Scores
  factor.scores <- ysvd$u
  fsncol <- ncol(factor.scores)
  for (j in 1:fsncol){
    factor.scores[,j] <- ysvd$u[,j] * -sqrt.ysvd[j]
  }
  #Add the column names
  PC.factorScore <- list(NULL, paste("Factor.scores.", 1:fsncol, sep=""))
  dimnames(factor.scores) <- PC.factorScore  
  #Add the row names
  fsnrow <- nrow(factor.scores)
  rownames(factor.scores) <- rownames(factor.scores, do.NULL=FALSE, prefix="Z.")
  
  #We multiply the fist column of the Factor scores with the first row of
  #Feature loading. This results in the predicted values
  FS <- subset(factor.scores, select = Factor.scores.1)
  #Create vector with the Loading values
  FL <- subset(feature.loading, select = Feature.Loading.1)
  TFL <- t(FL)
  FSTFL <- FS %*% TFL
  
  #Now we calculate the residuals
  Residuals <- x - FSTFL
  
  #Plot the vector coordinates
  plot <- plot(ysvd$v)
  
  ds = sum(x * x) #Data scatter
  mu1 <- round(ysvd$d[1] ^2) #Mu
  contr <- (mu1 / ds)*100 #Contribution of first component
  
  answer <- list(ysvd$d, sqr.ysvd, percent.ysvd, sqrt.ysvd ,ysvd$u, ysvd$v, ds, mu1, contr, x, transp.feature.loading, factor.scores, plot, TFL, FS, FSTFL, Residuals)
  names(answer)[[1]] <- "Singular.value"
  names(answer)[[2]] <- "Singular.value.squared"
  names(answer)[[3]] <- "Singular.value.percent"
  names(answer)[[4]] <- "Principal.component(Singular.value.sqrt)"
  names(answer)[[5]] <- "Scatter.coordinates"
  names(answer)[[6]] <- "PC.vector.coordinates"
  names(answer)[[7]] <- "Data.scatter"
  names(answer)[[8]] <- "Mu"
  names(answer)[[9]] <- "First.component.Contribution.%"
  names(answer)[[10]] <- "Original.dataset"
  names(answer)[[11]] <- "Feature.loads"
  names(answer)[[12]] <- "Factor.scores"
  names(answer)[[13]] <- "Plot"
  names(answer)[[14]] <- "PC1.Feature.loads"
  names(answer)[[15]] <- "PC1.Factor.Scores"
  names(answer)[[16]] <- "Predicted.values"
  names(answer)[[17]] <- "Residuals"
  return(answer)
}