source("/home/uu85g9/factor-analysis/mcem_fa_algorithm.R")
Y = readRDS("/data/users/uu85g9/factor_analysis-oldold/consolodation_may-24/processed_fa_scans.rds")
save_dir = "/data/users/uu85g9/factor-analysis/noms/mri/dpe_random_pxsigma/"
dir.create(save_dir, showWarnings=FALSE)
dir.create(paste0(save_dir, "intermediate_results/"), showWarnings=FALSE)
set.seed(1235)
lambda = matrix(runif(ncol(Y)*100, 0, 0.1), ncol(Y), 100)
v0s = c(0.1, 0.01, 0.001, 0.0001)
res = dpe_func(v0s=v0s, Y=Y, K=100, v1=10, do_stick_breaking=TRUE, tol=0.01, conv_param=TRUE, maxiter=2000,
               save_dir=paste0(save_dir, "intermediate_results/"), verbose=0, lambda=lambda,
               px_rotate=FALSE, varimax_every=-1, px_sigma=TRUE, nburn=200, n_mcmc=50, save_every=10)
saveRDS(res, paste0(save_dir, "res.rds"))
