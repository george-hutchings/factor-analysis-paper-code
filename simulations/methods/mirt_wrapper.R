# code to run veronika's method
# source("~/OneDrive/Documents/academic/24-25/factor-analysis/simulations/troubleshooting/compare-methods/veronika_method/probit_em_mirt.R")
# source("~/OneDrive/Documents/academic/24-25/factor-analysis/simulations/troubleshooting/compare-methods/veronika_method/probit_em_util.R")

source("~/factor-analysis/simulations/methods/veronika-files/probit_em_mirt.R")
source("~/factor-analysis/simulations/methods/veronika-files/probit_em_util.R")


mirt_wrapper = function(Y, K, ...){
  # taken from motivating example they used
  Y = Y-1
  stopifnot(sort(unique(as.numeric(Y)))==c(0,1))
# Initialize Loadings
large_k <- K
nitems <- ncol(Y)
loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
                             "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
lambda0_path <- c(0.2, 1, 5, 10, 20, 30, 40, 50)
lambda1 <- 0.2

# PX-EM
start = Sys.time()
px_em <- dynamic_posterior_exploration(data = Y, k = large_k, ibp_alpha = 2, mc_samples =50,
                                       ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                       max_iterations= 100, epsilon = 0.07, PX = TRUE, varimax = FALSE,
                                       loading_constraints= NULL, start = loading_starts_large,
                                       plot=FALSE, stop_rotation=100, random_state = 1, cores=1)

w = px_em$lambda0_40$gammas
lambda=px_em$lambda0_40$alphas
sf = diag(lambda%*%t(lambda) ) + 1
lambda_corr = lambda / sqrt(sf) 

return(list(list(w=w, lambda=lambda, lambda_corr=lambda_corr)))
}

