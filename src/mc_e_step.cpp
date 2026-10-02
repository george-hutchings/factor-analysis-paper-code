// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "utils.h"
using namespace arma;


//function(lambda, Z, alpha, group_indices, EZdTZd = colSums(Z**2), nsteps=50){
// [[Rcpp::export]]
Rcpp::List mce_step_rcpp( const mat &Y, const mat &lambda, mat &Z0, const vec &alpha, const Rcpp::List &group_indices, const vec &EZdTZd0, const uvec &binary, unsigned int nsteps=50)
{

  /* Useful variables that do not change each MCMC iteration */
  unsigned int K = lambda.n_cols;
  mat sigmaeta_inv = lambda.t() * lambda + diagmat(ones(K)); // A = t(lambda) %*% lambda + diag(1, K)
  mat cholsigma_eta = inv(trimatu( chol(sigmaeta_inv, "upper") ));
  mat sigmaeta = cholsigma_eta * cholsigma_eta.t();
  mat sigmaeta_lambda = sigmaeta * lambda.t();  // Matches R's sigmaeta_lambda = tcrossprod(sigmaeta, lambda)





  mat Z = Z0;
  unsigned int N = Z.n_rows;
  unsigned int D = Z.n_cols;


  mat eta(N,K,fill::zeros);
  mat EetaTeta(K,K,fill::zeros); /* E[etaTeta] */
  vec EZdTZd(D,fill::zeros);     /* E[z_d^TZ_d] GEORGE, not sure why set to uvec */
  mat EetaTZ(K,D,fill::zeros);   /* E[z_d^TZ_d] */
  mat Eeta(N,K,fill::zeros);     /* E[eta] */
  mat EZ(N,D,fill::zeros);       /* E[Z] */


  double logprior = 0.0;
  double reconstruction_error = 0.0;

  for (unsigned int b=0 ; b<nsteps ; b++ ){

    /* Sample eta */
    //samplescores_rcpp( eta, lambda, Z, alpha, Qchol );
    //samplescores_original_rcpp( eta, lambda, Z, alpha, sigmaeta, sigmaeta_lambda);
    samplescores_fixed_rcpp( eta, lambda, Z, alpha, cholsigma_eta, sigmaeta_lambda);
    EetaTeta += eta.t() * eta;
    Eeta += eta;

    /* Sample Z */
    logprior += sampleZ_rcpp( Y, lambda, eta, Z, group_indices, alpha, binary );
    logprior += eval_prior_eta( eta, lambda, alpha, Z );
    reconstruction_error += eval_reconstruction_error( eta, lambda, alpha, Z );

    /* Update running sums */
    EZdTZd += sum( square(Z), 0 ).t();
    EZ += Z;
    EetaTZ += eta.t() * Z;

  }

  /* Average over MCMC iterations */
  logprior /= nsteps;
  reconstruction_error /= nsteps;
  EetaTeta /= nsteps;
  Eeta /= nsteps;
  EZdTZd /= nsteps;
  EetaTZ /= nsteps;
  EZ /= nsteps;

  /* Output */
  return Rcpp::List::create(
    Rcpp::Named("Z") = Z,
    Rcpp::Named("Eeta") = Eeta,
    Rcpp::Named("EetaTeta") = EetaTeta,
    Rcpp::Named("EZ") = EZ,
    Rcpp::Named("EetaTZ") = EetaTZ,
    Rcpp::Named("EZdTZd") = EZdTZd,
    Rcpp::Named("logprior") = logprior,
    Rcpp::Named("reconstruction_error") = reconstruction_error
  );
}
