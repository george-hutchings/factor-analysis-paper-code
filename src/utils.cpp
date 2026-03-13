
#include "RcppArmadillo.h"
#include "rng.h"
#include "rtn1.h"
using namespace arma;


vec make_min_vals( const Rcpp::NumericVector &x )
{
  unsigned int n = x.size();
  vec vals(n+1);
  vals(0) = -INFINITY;
  for (unsigned int i=0 ; i<n ; i++ ) {
    vals(i+1) = x[i];
  }
  return(vals);
}
vec make_max_vals( const Rcpp::NumericVector &x )
{
  unsigned int n = x.size();
  vec vals(n);
  for (unsigned int i=1 ; i<n ; i++ ) {
    vals(i-1) = x[i];
  }
  vals(n-1) = INFINITY;
  return(vals);
}



// [[Rcpp::export]]
Rcpp::List range_by_group( const Rcpp::List &group_indices, const mat &Z, int col ) {
  
  int n = group_indices.size();
  Rcpp::NumericVector min_vals(n);
  Rcpp::NumericVector max_vals(n);
  
  for (int i = 0; i < n; i++) {
    Rcpp::IntegerVector idx = group_indices[i];
    double min_val = Z( idx[0] - 1, col ); // R is 1-based
    double max_val = min_val;
    
    for (int j = 1; j < idx.size(); j++) {
      double val = Z( idx[j] - 1, col );
      if (val < min_val) min_val = val;
      if (val > max_val) max_val = val;
    }
    
    min_vals[i] = min_val;
    max_vals[i] = max_val;
  }
  
  
  return Rcpp::List::create(Rcpp::Named("min_vals") = min_vals,
                            Rcpp::Named("max_vals") = max_vals);
}



// [[Rcpp::export]]
void samplescores_fixed_rcpp(mat &eta, const mat &lambda, const mat &Z, const vec &alpha, 
                             const mat &cholsigma_eta, const mat &sigmaeta_lambda) {
  unsigned int N = Z.n_rows;
  unsigned int K = lambda.n_cols;

  // Compute mu = sigmaeta_lambda * (Z - alpha)^T
  mat Z_centered = Z.each_row() - alpha.t();  // Subtract alpha from each column of Z
  mat mu = sigmaeta_lambda * Z_centered.t();

  // Generate noise: cholsigma_eta * N(0, I)
  mat noise(K, N, fill::randn);  // K x N matrix of standard normal noise
  mat tmp = cholsigma_eta * noise + mu;  // Transform noise and add mean

  // Assign transposed result to eta
  eta = tmp.t();  // Transpose to N x K
}

// [[Rcpp::export]]
void samplescores_rcpp( mat &eta, const mat &lambda, const mat &Z, const vec &alpha, const mat &Qchol )
{
  unsigned int N = Z.n_rows;
  unsigned int K = lambda.n_cols;
  vec tmp(K);
  
  for ( unsigned int n=0 ; n<N ; n++ ){
    
    tmp = lambda.t() * (Z.row(n).t()-alpha) ;
    eta.row(n) = rueMVnorm( tmp, Qchol ).t();
  }
}



// [[Rcpp::export]]
void samplescores_original_rcpp( mat &eta, const mat &lambda, const mat &Z, const vec &alpha, const mat &sigmaeta, const mat &sigmaeta_lambda)
{
  unsigned int K = lambda.n_cols;
  unsigned int N = eta.n_rows;
  
  mat noise(N,K,fill::randn);
  eta = noise * sigmaeta.t();
  
  for ( unsigned int n=0 ; n<N ; n++ ){
    eta.row(n) += (Z.row(n) - alpha.t()) * sigmaeta_lambda;
  }
  
}



// [[Rcpp::export]]
double sampleZ_rcpp( const mat &Y, const mat &lambda, const mat &eta, mat &Z, const Rcpp::List &group_indices, const vec &alpha, const uvec &binary )
{
  double logprior = 0.0;
  unsigned int n;
  unsigned int N = Z.n_rows;
  unsigned int D = Z.n_cols;
  
  mat a(N,D,fill::zeros);
  mat b=a, mu=eta*lambda.t();
  Rcpp::List tmp_list;
  vec min_vals, max_vals;
  uvec missing(N);

  
  for ( unsigned int d=0 ; d<D ; d++ ) {
    
    /* For each column, find min and max vals */ 
    tmp_list = range_by_group( group_indices[d], Z, d );
    min_vals = make_min_vals( tmp_list["max_vals"] );
    max_vals = make_max_vals( tmp_list["min_vals"] );
    
    /* Find indices of missing values */ 
    missing.fill(0);
    missing.elem( find_nonfinite(Y.col(d)) ).ones();
    
    /* Evaluate mu */
    mu.col(d) += alpha(d); 
    
    /* Repeat for each observation */ 
    for ( n=0 ; n<N ; n++) {
      
      if ( missing(n)==1 ) {
        
        /* First case: missing value */
        a(n,d) = -INFINITY;
        b(n,d) = INFINITY;
        
      } else {
        /* Second case: no missing value */
        if ( binary(d)==1 ) {
          /* Binary variable */
          if ( Y(n,d)==1 ) {
            a(n,d) = -INFINITY;
            b(n,d) = 0;
          } else {
            a(n,d) = 0;
            b(n,d) = INFINITY;
          }
        } else {
          /* Non-binary variable */
          a(n,d) = min_vals( Y(n,d)-1 );
          b(n,d) = max_vals( Y(n,d)-1 );
        }
      }
      
      /* Generate the z-value */
      Z.at(n,d) = rtn1( mu.at(n,d), 1.0, a.at(n,d), b.at(n,d) );

      /* Interval probability removed - belongs to sampler, not model likelihood */


    } /* End of loop for observations */
    
  } /* End of loop for variables */
  
  return(logprior);
}



double eval_prior_eta( const mat &eta, const mat &lambda, const vec &alpha, const mat &Z )
{
  unsigned int N = eta.n_rows;
  unsigned int D = Z.n_cols;
  unsigned int K = eta.n_cols;
  double logprior0 = 0.0;

  mat tmp = eta*lambda.t();
  unsigned int n, d;
  for ( n=0 ; n<N ; n++ ) {

    tmp.row(n) += alpha.t();

    /* Contribution from Z */
    for ( d=0 ; d<D ; d++ ) {
      logprior0 += R::dnorm( Z(n,d), tmp(n,d), 1.0, true );
    }

    /* Contribution from eta */
    for ( d=0 ; d<K ; d++ ) {
      logprior0 += R::dnorm( eta(n,d), 0.0, 1.0, true );
    }
  }

  return(logprior0);

}


double eval_reconstruction_error( const mat &eta, const mat &lambda, const vec &alpha, const mat &Z )
{
  unsigned int N = eta.n_rows;
  unsigned int D = Z.n_cols;
  double logprior0 = 0.0;

  mat tmp = eta*lambda.t();
  unsigned int n, d;
  for ( n=0 ; n<N ; n++ ) {

    tmp.row(n) += alpha.t();

    /* Contribution from Z only - no eta prior */
    for ( d=0 ; d<D ; d++ ) {
      logprior0 += R::dnorm( Z(n,d), tmp(n,d), 1.0, true );
    }
  }

  return(logprior0);

}
