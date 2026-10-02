
#include <RcppArmadillo.h>

using namespace arma;


vec rueMVnorm( const vec &b, const mat &Q ) {
  
  /* Remember that Q here is the lower Cholesky */ 
  /* Normally, the Cholesky is done within the function */ 
  /* Changed due to George's setup */ 
  
  int N = Q.n_cols;
  arma::vec b_new(N,fill::zeros);
  int i, j;
  
  
  /* Find v: Lv=b */
  arma::vec v(N); v.fill(0.0);
  v(0) = b(0)/Q(0,0);
  for (i=1 ; i<(N) ; i++) {
    v(i) = b(i); 
    for (j=0 ; j<i ; j++) {
      v(i) -= Q(i,j)*v(j);
    }
    v(i) /= Q(i,i);
  }
  
  
  /* Find m: Lm = v */
  arma::vec m(N); m.fill(0.0);
  m(N-1) = v(N-1)/Q(N-1,N-1);
  for (i=(N-2) ; i>=0 ; --i) {
    m(i) = v(i);
    for (j=(i+1) ; j<N ; j++) {
      m(i) -= Q(j,i)*m(j);
    }
    m(i) /= Q(i,i);
  }
  
  
  /* Generate z from N(0,I) */
  arma::vec z(N); z.fill(0);
  for (i=0 ; i<N ; i++) { 
    z(i) = R::rnorm(0.0,1.0);
  }
  
  
  /* Find w: Lw=z */
  arma::vec w(N); w.fill(0.0);
  w(N-1) = z(N-1)/Q(N-1,N-1);
  for (i=(N-2) ; i>=0 ; --i) {
    w(i) = z(i);
    for (j=(i+1) ; j<N ; j++) {
      w(i) -= Q(j,i)*w(j);
    }
    w(i) /= Q(i,i);
  }
  
  
  /* Add m and w */
  for (i=0 ; i<N ; i++) {
    b_new(i) = m(i) + w(i);
  }

  /* Output */
  return(b_new);
  
}