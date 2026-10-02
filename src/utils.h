

using namespace arma;


void samplescores_rcpp( mat &, const mat &, const mat &, const vec &, const mat & );

void samplescores_original_rcpp( mat &, const mat &, const mat &, const vec &, const mat &, const mat &);

void samplescores_fixed_rcpp( mat &, const mat &, const mat &, const vec &, const mat &, const mat &);

double sampleZ_rcpp( const mat &, const mat &, const mat &, mat &, const Rcpp::List &, const vec &, const uvec & );

double eval_prior_eta( const mat &, const mat &, const vec &, const mat & );

double eval_reconstruction_error( const mat &, const mat &, const vec &, const mat & );
