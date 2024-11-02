#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2018  Wenjie Wang <wjwang.stat@gmail.com>
//
// This file is part of the R package intsurv.
//
// The R package intsurv is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package intsurv is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
 
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>
using namespace arma;
using namespace Rcpp;
using namespace std;
 
// cumulative sum in possibly reverse order
// [[Rcpp::export]]
inline arma::vec cum_sum(const arma::vec& x,
                         const bool reversely = false)
{
  // if cumsum reversely
  if (reversely) {
    const unsigned long int n_x {x.n_rows};
    arma::vec res {arma::zeros(n_x)};
    double tmp {0.0};
    for (size_t i {1}; i <= n_x; ++i) {
      tmp += x[n_x - i];
      res[n_x - i] = tmp;
    }
    return res;
  }
  // otherwise, using arma::cumsum
  return arma::cumsum(x);
}


// column-wise cumulative sum in possibly reverse order
// [[Rcpp::export]]
inline arma::mat cum_sum_cols(const arma::mat& x,
                              const bool reversely = false)
{
  // if cumsum reversely
  if (reversely) {
    const unsigned long int n_x = x.n_rows;
    arma::mat tmp {arma::zeros(1, x.n_cols)};
    arma::mat res {x};
    for (size_t i {1}; i <= n_x; ++i) {
      tmp += x.row(n_x - i);
      res.row(n_x - i) = tmp;
    }
    return res;
  }
  // otherwise, using arma::cumsum
  return arma::cumsum(x, 0);
}







// [[Rcpp::export]]
arma::mat Xotimes2(arma::mat X){
  // t(apply(X, 1, function(x) c(x %o% x )))
  int px = X.n_cols, n = X.n_rows, i;   
  arma::mat X2(n, px*px, fill::zeros), m; 
  arma::rowvec Xi; 
  
  for(i=0; i<n; i++){
    Xi = X.row(i);
    m = Xi.t() * Xi;
    X2.row(i) = m.as_row(); 
  }
  
  return(X2);
}
