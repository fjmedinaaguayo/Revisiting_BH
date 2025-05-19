// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double log_factorial_cpp(int n) {
  return lgamma(n+1);
}

NumericVector my_dnorm( NumericVector x, NumericVector means, NumericVector sds, bool isLog){
  int d = x.size() ;
  NumericVector res(d) ;
  for( int i=0; i<d; i++) res[i] = R::dnorm( x[i], means[i], sds[i], isLog) ;
  return res ;
}

// [[Rcpp::export]]
double logQ_L_cpp(NumericVector x,NumericVector mu_l,NumericVector sd_prop, bool ordrd = false) {

  int d = x.size() ;
  int k = d/2;
  double out=sum(my_dnorm(x,mu_l,sd_prop,true));
  
  if(ordrd==true)
    out = out + log_factorial_cpp(k);
  
  return out;
}

// [[Rcpp::export]]
NumericVector logQ_L_vec_cpp(NumericVector x, NumericMatrix mu_vec_sub, NumericMatrix sigma_vec_sub, bool ordrd = false){
  
  int n=mu_vec_sub.rows();
  NumericVector out(n);
  
  for(int i=0; i<n; i++){
    
    NumericVector mu_l=mu_vec_sub(i,_);
    NumericVector sigma_l=sigma_vec_sub(i,_);
    out[i]=logQ_L_cpp(x,mu_l,sigma_l,ordrd);
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector logGamma_den_cpp(NumericVector x, NumericMatrix mu_vec_sub, NumericVector N_vec, NumericVector sd_prop, double beta_t){
  
  int n=mu_vec_sub.rows();
  NumericVector out(n);
  
  for(int i=0; i<n; i++){
    
    NumericVector mu_l=mu_vec_sub(i,_);
    out[i]=beta_t*logQ_L_cpp(x,mu_l,sd_prop) + log(N_vec[i]);
  }
  
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
