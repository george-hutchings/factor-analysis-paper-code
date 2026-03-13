#change initialisation, e step in one function, change expectations needed, and change in report change main so that is working

# CHECKPOINT RESUME FUNCTIONALITY:
# - mcem_algorithm: Pass resume_from="path/to/iter100.rds" to continue from iteration 101
# - dpe_func: Pass resume_dpe=TRUE to auto-detect completed DPE stages and resume between stages
#   Note: DPE resume checks for completed stages (finds any checkpoint = stage done)
#   and starts from first incomplete stage. Does not continue partial stages to higher maxiter.
#   Always provide same v0s vector when resuming (v0 values checked with warning if mismatch).

library(RcppTN)
library(nloptr)

#for inverse gamma operations
library(invgamma)

# for faster calculation of boundaries
library(Rcpp)

# C++ stuff
library(Rcpp)
library(RcppArmadillo)


# Use Q-function corrected version
sourceCpp(file.path(dirname(sys.frame(1)$ofile), "src/mc_e_step.cpp"))

mcem_algorithm = function(Y, K, v0=1e-4, v1=10, tol=0.05, maxiter=2000, save_dir=NULL, verbose=FALSE, lambda=NULL, Z=NULL, alpha=NULL, thetas=NULL, conv_param = FALSE, do_stick_breaking=TRUE,nburn=100, n_mcmc=50, px_sigma=TRUE, px_rotate=TRUE, varimax_every=30, resume_from=NULL, save_every=10, ibp_alpha=2){
  # conv_param determines whether to use the parameter differences as convergence condition
  
  missing_mask = is.na(Y)
  
  # useful variables
  N = nrow(Y)
  D = ncol(Y)
  
  # vector indicating whether dimensions are binary or not
  binary = rep(FALSE, D)
  # function to map classes to groups 1,2,... which is useful for indexing
  for (d in seq_len(D)){
    uniq_Y = sort(unique(Y[,d])) # ignores NA
    if(length(uniq_Y)==2){
      binary[d]=TRUE
    } 
    Y[,d] = match(Y[,d], uniq_Y) # ignores NA
  }
  binary_rcpp = 1*binary
  # produces a list of form x[[dim]] = list s.t list[[i]] = row indicies of Y[,dim] that are in group i
  # this is useful for the extended rank likelihood
  group_indices <- lapply(seq_len(D), function(d) split(seq_len(nrow(Y)), Y[, d])) # ignores NA
  
  # Resume from checkpoint if provided
  start_iter = 1
  prev_parameter_diff = prev_q_function = NULL
  if (!is.null(resume_from)) {
    if (!file.exists(resume_from)) stop(paste("Resume file not found:", resume_from))
    checkpoint = readRDS(resume_from)
    if (!all(c("lambda", "Z", "alpha", "thetas", "iterations") %in% names(checkpoint))) {
      stop("Invalid checkpoint file")
    }
    lambda = checkpoint$lambda
    Z = checkpoint$Z
    alpha = checkpoint$alpha
    thetas = checkpoint$thetas
    start_iter = checkpoint$iterations + 1
    prev_parameter_diff = checkpoint$parameter_diff
    prev_q_function = checkpoint$q_function
    print(paste("Resuming from iteration", checkpoint$iterations))
    nburn = 0
  }
  
  # useful cpp function that makes compute much faster(~x20) when computing the boundaries (when there are lots of them)
  cppFunction('
Rcpp::List range_by_group(Rcpp::List group_indices, Rcpp::NumericMatrix Z, int col) {
    int n = group_indices.size();
    Rcpp::NumericVector min_vals(n);
    Rcpp::NumericVector max_vals(n);
    
    for (int i = 0; i < n; i++) {
        Rcpp::IntegerVector idx = group_indices[i];
        double min_val = Z(idx[0] - 1, col); // R is 1-based
        double max_val = min_val;

        for (int j = 1; j < idx.size(); j++) {
            double val = Z(idx[j] - 1, col);
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }

        min_vals[i] = min_val;
        max_vals[i] = max_val;
    }
    
    return Rcpp::List::create(Rcpp::Named("min_vals") = min_vals,
                              Rcpp::Named("max_vals") = max_vals);
}')
  
  
  # Initialise values
  
  # initial Z
  rank_with_na = function(x){
    # finds ranks, assigning na random rank
    tmp = rank(x, na.last='keep')
    tmp[is.na(tmp)] = sample(tmp[!is.na(tmp)], sum(is.na(tmp)), replace=TRUE)
    tmp = rank(tmp, ties.method="average")
    return(tmp)
  }
  
  if (is.null(Z)){
    tmp = apply(Y, 2, rank_with_na)/(N+1) # ranks
    Z = qnorm(tmp)
    Z = scale(Z)
  }
  
  
  # initialising alpha, if binary initialise from base rate (just initialisation)
  if (is.null(alpha)){
    alpha = rep.int(0L, D)
    alpha[binary] = qnorm(colMeans(Y[,binary, drop=FALSE], na.rm=TRUE) - 1)
  }
  ## initialising lambda - this is done using PCA 
  if (is.null(lambda)){
    # 1. initialise lambda using pca
    lambda = matrix(0, D, K) 
    tmp =  prcomp(Z, rank.=K, center = TRUE,scale. = TRUE) #pca components limited to D, K
    lambda[, 1:min(D,K)] = t( t(tmp$rotation)*tmp$sdev[1:min(D,K)] )
    
    # 2. initialise any excess columns randomly so that they are smaller than the principal components
    # rnorm(D*K, sd=1/sqrt(D) gives columns of norm = 1, then we scale them by the smallest e/val *0.5
    if(D<K){
      lambda[, D:K] = rnorm(D*(K-D+1), sd=1/sqrt(D))*tmp$sdev[min(D,K)]
    }
    
  }else{
    # if lambda is given checking the the dimensions align with K an D (from Y)
    stopifnot(nrow(lambda)==D)
    stopifnot(ncol(lambda)==K)
  }
  
  
  # initialise w (inclusion probabilities)
  w = matrix(0.5, D, K)
  # initialise thetas
  if (is.null(thetas)){
    #thetas = seq(1, 0, length.out=K+2)[2:(K+1)] 
    thetas<-sort(rbeta(K,1,1),decreasing=TRUE)
    
  }
  
  # calculate here
  # vector with elements E[Z_d^TZ_d]
  EZdTZd = colSums(Z**2)
  
  
  # calculates the relevant expectations for w ( E- step )
  update_w = function(lambda, thetas, v1, v0){
    tmp =  t(dnorm( lambda, sd=sqrt(v1) ))*thetas # in KxD
    t(
      tmp / (tmp + t(dnorm( lambda, sd=sqrt(v0) ))*(1 - thetas) )
    )
  }
  
  update_alpha = function(lambda, EZ, Eeta, N=nrow(Eeta)){
    #E[ (Zd - eta lambdad)^T ONES ]
    EZminusetaLambda = colSums( EZ - tcrossprod(Eeta, lambda) )
    EZminusetaLambda/(N + 1)
  }
  
  # cpp update lambda function
  cppFunction(depends = "RcppArmadillo", code = '
arma::mat update_lambda_cpp(const arma::mat& EetaTZ, 
                            const arma::mat& Eeta,
                            const arma::vec& alpha,
                            const arma::mat& w,
                            double v0,
                            double v1,
                            const arma::mat& EetaTeta) {
    
    unsigned int D = w.n_rows;
    
    // Calculate E[etaT(Z-alpha)]
    arma::mat EetaTZminusalpha = EetaTZ - sum(Eeta, 0).t() * alpha.t();
    
    // Calculate EHinv
    arma::mat EHinv = trans( (w * (v0 - v1) + v1) / (v1 * v0) );
    
    // Initialize output matrix
    arma::mat lambda(D, w.n_cols);
    
    // preallocate sigmad where we just update diagonal each iteration
    arma::mat sigmad = EetaTeta;
    arma::vec diag_sigmad = sigmad.diag();
    
    // For each d, calculate lambda[d,]
    for (unsigned int d = 0; d < D; d++) {
        // Create sigmad = EetaTeta + diag(EHinv.col(d))
        sigmad.diag() = diag_sigmad + EHinv.col(d);
        
        // Use solve() instead of inv_sympd() for better performance when multiplying immediately
        // This solves the system sigmad * x = EetaTZminusalpha.col(d) directly
        // this seems to provide massive speed up (esp force_sympd)
        lambda.row(d) = trans(solve(sigmad, EetaTZminusalpha.col(d), arma::solve_opts::likely_sympd));
        
    }
    
    return lambda;
}
')
  
  # prior parameters for sigma ( small value = uninformative )
  sigma_alpha0 = sigma_beta0 = 0.001
  # posterior alpha parameter
  sigma_alphaplus1 = sigma_alpha0 + N/2 + 1
  update_sigma = function(EZdTZd, alpha, lambda, EetaTeta, EZ, Eeta, sigma_beta0, sigma_alphaplus1, estimate_sigma, N = nrow(EZ)){
    #E[ (Zd - alphad - eta lambdad)^T(Zd - alphad - eta lambdad) ]
    EZminusalphaminsetalambda2 = EZdTZd + N*(alpha**2) + rowSums( ( lambda%*%EetaTeta )*lambda) -
      2*colSums(EZ)*alpha + 2*rowSums(tcrossprod(lambda, Eeta))*alpha - 2*colSums(t(lambda)*EetaTZ)
    b = sigma_beta0 + EZminusalphaminsetalambda2*0.5
    sigma = b/sigma_alphaplus1
    sigma
  }
  
  ## Functions for non-linear optimiser for theta.
  # NB these functions are undefined for x=0, x=1 so 1e-9 is included to avoid
  # numerical instability
  
  # objective function
  theta_Q = function(x, sum_w, D, ibp_alpha){
    if (any(is.na(x)) ||  x[1] == 1){
      return( NaN )
    }else{
      return( -sum( sum_w*log(x) + log(1-x)*(D - sum_w) ) - (ibp_alpha-1)*log(x[length(x)]) ) 
    }
  }
  
  theta_grad_Q = function(x, sum_w, D, ibp_alpha){
    ans = -( sum_w/(x) - (D - sum_w)/(1-x) )
    ans[length(x)] = (ans[length(x)] + (1-ibp_alpha)/x[length(x)])
    ans
  }
  theta_constraint_Q = function(x, sum_w, D, ibp_alpha){
    x[-1]-x[-length(x)] 
  }
  theta_jac_constraint_Q = function(x, sum_w, D, ibp_alpha){
    lenx = length(x)
    mat = diag(-1, nrow=lenx-1,ncol=lenx) #k-1xk
    tmp = seq.int(length.out=lenx-1)
    mat[cbind(tmp,tmp+1)] = 1
    mat
  }
  
  if (do_stick_breaking){
    update_theta = function(w, thetas, ibp_alpha, theta_Q, theta_grad_Q, theta_constraint_Q, theta_jac_constraint_Q, D = nrow(w)){
      sum_w = colSums(w)
      #theta0 = (thetas*1e6 +0.5) / (1e6 + 1) # weighted mean to ensure x0 is valid
      theta0<-sort(rbeta(K,1,1),decreasing=TRUE)
      thetas = nloptr(x0=theta0,eval_f=theta_Q, eval_grad_f=theta_grad_Q, eval_g_ineq=theta_constraint_Q,
                      eval_jac_g_ineq=theta_jac_constraint_Q,
                      opts=list("algorithm"="NLOPT_LD_MMA","check_derivatives"=F, "xtol_rel"=10^-10),
                      lb=rep_len(0,K),ub=rep_len(1,K),
                      D=D,
                      sum_w=sum_w,
                      ibp_alpha=ibp_alpha)$solution
      thetas
    }
  }else{
    update_theta = function(w, thetas, ibp_alpha, theta_Q, theta_grad_Q, theta_constraint_Q, theta_jac_constraint_Q, D = nrow(w)){
      colMeans(w)
    }
  }
  
  # burn in
  Z = mce_step_rcpp( Y=Y, lambda=lambda, Z0=Z, alpha=alpha, group_indices, EZdTZd0=EZdTZd, binary_rcpp, nsteps=nburn)$Z
  
  
  # Initialize or extend convergence arrays
  # Ereconstruction = E[log N(Z | Lambda*eta + alpha, I)] = -0.5 * E[||Z - Lambda*eta - alpha||^2] + const
  if (!is.null(resume_from)) {
    parameter_diff = q_function = q_monotonic = Ereconstruction = rep(NA, maxiter)
    parameter_diff[1:length(prev_parameter_diff)] = prev_parameter_diff
    q_function[1:length(prev_q_function)] = prev_q_function
    if (!is.null(checkpoint$q_monotonic)) {
      q_monotonic[1:length(checkpoint$q_monotonic)] = checkpoint$q_monotonic
    }
    lambda_old = lambda  # Use checkpoint lambda to avoid Inf diff
  } else {
    parameter_diff = q_function = q_monotonic = Ereconstruction = rep(NA, maxiter)
    lambda_old = matrix(Inf, D, K)
  }

  # Save initial state (iter=0) if save_dir provided
  if(!is.null(save_dir) && start_iter == 1){
    sf = 1/sqrt( rowSums(lambda**2) + 1 )
    saveRDS(
      list(alpha = alpha, lambda = lambda, lambda_corr=lambda*sf, w=w, thetas=thetas, sigma=NULL, iterations=0, v0=v0, v1=v1, parameter_diff=numeric(0), q_function=numeric(0), q_monotonic=numeric(0), Ereconstruction=numeric(0), Y=Y, Z=Z),
      file=paste0(save_dir, 'iter0.rds'))
  }

  start_time = Sys.time()
  for (iter in start_iter:maxiter){
    
    ## E Step
    mce_step_res = mce_step_rcpp( Y=Y, lambda=lambda, Z0=Z, alpha=alpha, group_indices, EZdTZd0=EZdTZd, binary_rcpp, nsteps=n_mcmc)
    EetaTeta = mce_step_res$EetaTeta
    Eeta = mce_step_res$Eeta
    
    Z = mce_step_res$Z
    EZ = mce_step_res$EZ
    
    
    EetaTZ = mce_step_res$EetaTZ

    EZdTZd = mce_step_res$EZdTZd

    Ereconstruction[iter] = mce_step_res$reconstruction_error

    # w
    w = update_w(lambda, thetas, v1, v0)

    # E[log p(Lambda | w)] - spike-and-slab prior (was missing in old version)
    log_pdf_v1 = dnorm(lambda, 0, sqrt(v1), log=TRUE)
    log_pdf_v0 = dnorm(lambda, 0, sqrt(v0), log=TRUE)
    expected_log_lambda = sum( w * log_pdf_v1 + (1 - w) * log_pdf_v0 )

    # Q(θ_old|θ_old): Expected complete-data log-posterior before M-step
    # NB: C++ logprior excludes interval probability (Gibbs sampler artifact, not model likelihood)
    q_current = mce_step_res$logprior +      # E[log p(Z|eta,Lambda,alpha) + log p(eta)]
      expected_log_lambda +                   # E[log p(Lambda|w)]
      sum(dnorm(alpha, log=TRUE)) +           # log p(alpha)
      sum(dbeta(pmax( pmin( thetas / c(1, thetas[-length(thetas)]), 1), 1e-9), ibp_alpha, 1, log = TRUE))  # log p(theta)
    # E[log p(w|theta)] - coefficient is D (one per loading), not 1
    for (k in seq_len(K)){
      q_current = q_current + sum( w[,k]*log( (thetas[k]+1e-9)/(1-thetas[k]+1e-9) ) ) + D * log(1 - thetas[k] + 1e-9)
    }
    q_function[iter] = q_current
    
    
    
    ## M-step
    lambda = update_lambda_cpp(EetaTZ, Eeta, alpha, w, v0, v1, EetaTeta)
    
    thetas = update_theta(w, thetas, ibp_alpha, theta_Q, theta_grad_Q, theta_constraint_Q, theta_jac_constraint_Q, D)
    
    alpha = update_alpha(lambda, EZ, Eeta, N)

    # ================================================================
    # Q(θ_new | θ_old) - GUARANTEED TO INCREASE BY EM THEORY
    # Uses: old expectations (EetaTeta, Eeta, EZ, EetaTZ, EZdTZd)
    # Uses: old w (posterior inclusion probabilities from E-step)
    # Evaluates at: new parameters (lambda, alpha, thetas)
    # ================================================================

    # Term 1: E[log p(Z | η, λ_new, α_new)] - reconstruction
    # E[||Z - ηλ^T - α||²] expanded using old expectations, new params
    E_ZZ = sum(EZdTZd)
    E_Z_eta_lambda = sum(lambda * t(EetaTZ))
    E_Z_alpha = sum(colSums(EZ) * alpha)
    E_eta_lambda_sq = sum(diag(crossprod(lambda) %*% EetaTeta))
    E_eta_lambda_alpha = sum(alpha * drop(lambda %*% colSums(Eeta)))
    E_alpha_sq = N * sum(alpha^2)

    E_resid_sq = E_ZZ - 2*E_Z_eta_lambda - 2*E_Z_alpha +
                 E_eta_lambda_sq + 2*E_eta_lambda_alpha + E_alpha_sq
    log_recon = -0.5 * E_resid_sq - 0.5 * N * D * log(2*pi)

    # Term 2: E[log p(η)] - eta prior (same for old/new, only uses expectations)
    log_eta = -0.5 * sum(diag(EetaTeta)) - 0.5 * N * K * log(2*pi)

    # Term 3: E[log p(λ_new | w)] - spike-and-slab at new lambda, old w
    log_pdf_v1_new = dnorm(lambda, 0, sqrt(v1), log=TRUE)
    log_pdf_v0_new = dnorm(lambda, 0, sqrt(v0), log=TRUE)
    log_lambda = sum(w * log_pdf_v1_new + (1 - w) * log_pdf_v0_new)

    # Term 4: log p(α_new) - alpha prior at new alpha
    log_alpha = sum(dnorm(alpha, 0, 1, log=TRUE))

    # Term 5: E[log p(w | θ_new)] - Bernoulli prior, old w, new theta
    log_w = 0
    for (k in seq_len(K)) {
      log_w = log_w + sum(w[,k] * log(thetas[k] + 1e-9)) +
                      sum((1 - w[,k]) * log(1 - thetas[k] + 1e-9))
    }

    # Term 6: log p(θ_new) - IBP stick-breaking prior at new theta
    nu = thetas / c(1, thetas[-length(thetas)])
    nu = pmax(pmin(nu, 1 - 1e-9), 1e-9)
    log_theta = sum(dbeta(nu, ibp_alpha, 1, log=TRUE))

    q_monotonic[iter] = log_recon + log_eta + log_lambda + log_alpha + log_w + log_theta

    if (px_sigma){
    sigma = update_sigma(EZdTZd, alpha, lambda, EetaTeta, EZ, Eeta, sigma_beta0, sigma_alphaplus1, estimate_sigma, N)
    sigma = drop(sigma)
    lambda = lambda/sqrt(sigma)
    alpha = alpha/sqrt(sigma)
    # Scale Z for Gibbs warm start (boundaries recalculated from Z in C++)
    Z = sweep(Z, 2, sqrt(sigma), FUN = "/")
    # Scale EZdTZd for consistency (not strictly needed - C++ recomputes from Z)
    EZdTZd = EZdTZd / sigma
    }

    if (px_rotate){
    a_mat <- mce_step_res$EetaTeta/N
    lambda <- tcrossprod(lambda, chol(a_mat))
    }

    if (varimax_every>0){
    if(iter %% varimax_every == 0){
      lambda=matrix(varimax(lambda+1e-6)$loadings,D,K)
    }
    }
    
    
    
    
    # convergence conditions
    parameter_diff[iter] = max( abs(abs(lambda) - abs(lambda_old)) )
    lambda_old=lambda
    
    if (verbose!=-1){
      print(paste('Iteration:', iter, 'Current parameter diff:', parameter_diff[iter]))
      print(paste('Iteration:', iter, 'Q(θ_old|θ_old):', round(q_function[iter], 2),
                  'Q(θ_new|θ_old):', round(q_monotonic[iter], 2),
                  'M-step gain:', round(q_monotonic[iter] - q_function[iter], 2)))
      print('Avg time per iteration:')
      print( difftime(Sys.time(), start_time, units="mins")/iter )
    }
    
    if (conv_param){
      if (iter > 1 && parameter_diff[iter] < tol){
        print(paste("Converged in", iter, "iterations"))
        break
      }
    } else {
      if (iter > 1) {
        if (!is.finite(q_function[iter]) || !is.finite(q_function[iter-1])) {
          print(paste("Infinite Q-function detected at iteration", iter))
          next  # Skip to next iteration
        } else if (abs(q_function[iter] - q_function[iter-1]) < tol) {
          print(paste("Converged in", iter, "iterations"))
          break
        }
      }
    }
    
    
    if (iter %% save_every == 0){
      if(!is.null(save_dir)){
        sf = 1/sqrt( rowSums(lambda**2) + 1 )
        sigma_save = if(px_sigma && exists("sigma")) sigma else NULL
        saveRDS(
          list(alpha = alpha, lambda = lambda, lambda_corr=lambda*sf, w=w, thetas=thetas, sigma=sigma_save, iterations=iter, v0=v0, v1=v1, parameter_diff=parameter_diff[1:iter], q_function=q_function[1:iter], q_monotonic=q_monotonic[1:iter], Ereconstruction=Ereconstruction[1:iter], Y=Y, Z=Z),
          file=paste0(save_dir, 'iter', iter, '.rds'))
      }
    }
    
    if (verbose==1){
      if(iter%%5 == 1){
        print(paste('Iteration:', iter))
        print('Lambda:')
        sf = 1/sqrt( rowSums(lambda**2) + 1 )
        print(lambda * sf)
        print('Theta:')
        print(thetas)
        print('W:')
        print(w)
        print('Alpha:')
        print(alpha)
        print(strrep('-', 80))
      }
    }
    
  }
  
  
  sf = 1/sqrt( rowSums(lambda**2) + 1 )
  sigma_save = if(px_sigma && exists("sigma")) sigma else NULL

  return(list(alpha = alpha, lambda = lambda, lambda_corr=lambda*sf, w=w, thetas=thetas, sigma=sigma_save, iterations=iter, v0=v0, v1=v1, parameter_diff=parameter_diff[1:iter], q_function=q_function[1:iter], q_monotonic=q_monotonic[1:iter], Ereconstruction=Ereconstruction[1:iter], Z=Z))
}


# Helper to find latest checkpoint for a DPE stage
find_latest_checkpoint = function(save_dir, dpe_stage) {
  if (is.null(save_dir)) return(NULL)
  dir_path = dirname(save_dir)
  if (dir_path == ".") dir_path = getwd()
  files = list.files(dir_path, pattern = paste0("dpe", dpe_stage, "_iter[0-9]+\\.rds$"), full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files[which.max(as.integer(sub(".*iter([0-9]+)\\.rds$", "\\1", files)))]
}

dpe_func = function(v0s, save_dir=NULL, lambda= NULL, nburn=150, resume_dpe=FALSE, ...){
  lambda_res = lambda
  Z_res = alpha_res = theta_res= NULL
  results = list()
  n_dpe = length(v0s)
  start_dpe = 1
  
  # Check for existing DPE results to resume from
  if (resume_dpe && !is.null(save_dir)) {
    for (i in seq_len(n_dpe)) {
      checkpoint_file = find_latest_checkpoint(save_dir, i)
      if (is.null(checkpoint_file)) {
        start_dpe = i
        if (i > 1 && !is.null(prev_checkpoint <- find_latest_checkpoint(save_dir, i-1))) {
          prev = readRDS(prev_checkpoint)
          if (!is.null(prev$v0) && prev$v0 != v0s[i-1]) {
            warning(paste("v0 mismatch: checkpoint has v0=", prev$v0, "but expected", v0s[i-1]))
          }
          lambda_res = matrix(varimax(prev$lambda+1e-6)$loadings, nrow(prev$lambda), ncol(prev$lambda))
          Z_res = prev$Z
          alpha_res = prev$alpha
          theta_res = prev$thetas
          nburn = 100
        }
        break
      }
      results[[i]] = readRDS(checkpoint_file)
      if (!is.null(results[[i]]$v0) && results[[i]]$v0 != v0s[i]) {
        warning(paste("v0 mismatch at stage", i, ": checkpoint has v0=", results[[i]]$v0, "but expected", v0s[i]))
      }
    }
    if (start_dpe > 1) print(paste("Resuming DPE from stage", start_dpe, "of", n_dpe))
  }
  
  print('Starting DPE')
  for (i in start_dpe:n_dpe){
    save_dir_dpei = if(!is.null(save_dir)) paste0(save_dir, "dpe", i, "_") else NULL
    resume_file = if(resume_dpe && !is.null(save_dir)) find_latest_checkpoint(save_dir, i) else NULL
    
    a = Sys.time()
    print(paste0('Starting DPE ', i, ' of ', n_dpe))
    tmp = mcem_algorithm(v0=v0s[i], lambda = lambda_res, Z = Z_res, alpha = alpha_res, thetas=theta_res, nburn=nburn, save_dir = save_dir_dpei, resume_from=resume_file, ... )
    results[[i]] = tmp
    D=nrow(tmp$lambda)
    K=ncol(tmp$lambda)
    lambda_res=matrix(varimax(tmp$lambda+1e-6)$loadings,D,K)
    Z_res = tmp$Z
    alpha_res = tmp$alpha
    theta_res = tmp$thetas
    nburn=100
    
    tmp = difftime( Sys.time(), a, units = 'mins')
    print(paste('DPE ', i, 'Complete in ', tmp, 'mins'))
  }
  return(results)
}
# lambda = t(t(lambda)*sign(colSums(lambda)))

