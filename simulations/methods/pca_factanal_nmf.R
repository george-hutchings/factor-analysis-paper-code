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
  # Scale component k by sdev[k], giving variable-component correlations. rank.=K
  # truncates rotation to D x K but leaves sdev at its full length min(n, p), so
  # `tmp$rotation * tmp$sdev` recycles column-major and scales row i by sdev[i]
  # instead. sweep over the column margin is explicit and cannot recycle.
  tmp = unclass(varimax(sweep(tmp$rotation, 2, tmp$sdev[seq_len(K)], "*"))$loadings)
  tmp = tmp[, order(colSums(abs(tmp)), decreasing=TRUE) ]
  tmp = t( t(tmp)*sign(colSums(tmp)))
  tmp # lambda corr
}

my_nmf = function(Y, K, ...){
  stopifnot(Y>0)
  
  # scaling so the data has unit variance and lambda corresponds to correlation
  Y = t( t(Y)/sqrt(diag(var(Y))) )
  
  # method='SNMF/L' puts the sparsity on W (the scores), so the loadings t(H) are not sparse
  tmp = nmf(Y, K, method='SNMF/L')
  error = Y - tmp@fit@W%*%tmp@fit@H

  lambda = t(tmp@fit@H)
  lambda = lambda[, order(colSums(abs(lambda)), decreasing=TRUE) ]
  lambda
}
