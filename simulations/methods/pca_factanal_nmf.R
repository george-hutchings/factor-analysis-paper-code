library(NMF)
#comparing methods
#factanal method
my_factanal = function(Y, K, ...){
  tmp = factanal(Y, K, nstart=20)
  # divide by uniquness to correspond to correlation
  lambda = unclass(tmp$loadings) #/sqrt(tmp$uniquenesses))
  lambda = lambda[, order(colSums(abs(lambda)), decreasing=TRUE) ]
  lambda = t( t(lambda)*sign(colSums(lambda)))
  lambda # lambda corr
}


#basic method PCA + varimax
pca_varimax = function(Y, K, ...){
  Y = scale(Y) # scaling so the data has unit variance
  tmp = prcomp(Y, rank.=K)
  tmp = unclass(varimax(tmp$rotation*tmp$sdev)$loadings)
  tmp = tmp[, order(colSums(abs(tmp)), decreasing=TRUE) ]
  tmp = t( t(tmp)*sign(colSums(tmp)))
  tmp # lambda corr
}

my_nmf = function(Y, K, ...){
  stopifnot(Y>0)
  
  # scaling so the data has unit variance and lambda corresponds to correlation
  Y = t( t(Y)/sqrt(diag(var(Y))) )
  
  # method='SNMF/L' is sparse loading matrix
  tmp = nmf(Y, K, method='SNMF/L')
  error = Y - tmp@fit@W%*%tmp@fit@H

  lambda = t(tmp@fit@H)
  lambda = lambda[, order(colSums(abs(lambda)), decreasing=TRUE) ]
  lambda
}
