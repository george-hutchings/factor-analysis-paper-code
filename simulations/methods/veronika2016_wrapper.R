# ***** Factor regression simulated example ****** #


source("~/factor-analysis/simulations/methods/FACTOR_CODE_update.R")

library(mvtnorm)
library(partitions)
library(nloptr)
library(glmnet)

veronika2016 = function(..., varimax_all = FALSE){
  # K, Y defined
  p=D = ncol(Y)
  alpha<-1/ncol(Y)
  lambda1<-0.001
  epsilon<-0.05 # default is 0.05
  startB = matrix(rnorm(D*K),D,K)
  
  # Each stage is warm-started from the previous one, as in R&G's EXAMPLE_FACTOR_ROTATE.R.
  # varimax_all rotates on every stage: TRUE in the high-dimensional study, where the
  # lambda0 = 20 stage does not converge without it, and FALSE in scenario 1, matching
  # our method, which is not rotated within stages.
  vmx = varimax_all

  print('start dpe 1')
  start<-list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))
  
  #Y,lambda0,lambda1,start,K,epsilon,alpha,PX,approximate,stop,varimax,plot=TRUE
  lambda0<-5
  result_5_2<-FACTOR_ROTATE(Y,lambda0,lambda1,start,K,epsilon,alpha,TRUE,FALSE,100,vmx)
  
  print('start dpe 2')
  lambda0<-10
  result_10_2<-FACTOR_ROTATE(Y,lambda0,lambda1,result_5_2,K,epsilon,alpha,TRUE,FALSE,100,vmx)
  
  print('start dpe 3')
  lambda0<-20
  result_20_2<-FACTOR_ROTATE(Y,lambda0,lambda1,result_10_2,K,epsilon,alpha,TRUE,FALSE,100,vmx)
  
  print('start dpe 4')
  lambda0<-30
  result_30_2<-FACTOR_ROTATE(Y,lambda0,lambda1,result_20_2,K,epsilon,alpha,TRUE,FALSE,100,vmx)
  
  lambda0<-40
  result_40_2<-FACTOR_ROTATE(Y,lambda0,lambda1,result_30_2,K,epsilon,alpha,TRUE,FALSE,100,vmx)
  
  result = result_40_2
  w = result$P_star
  lambda = result$B
  
  sf = 1/sqrt( diag(lambda%*%t(lambda) ) + result$sigma)
  lambda_corr = sf*lambda
  
  list(w=w, lambda=lambda, lambda_corr=lambda_corr)
}


