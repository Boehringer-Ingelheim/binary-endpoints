#include <Rcpp.h>
using namespace Rcpp;

bool check_limit(NumericVector lprob, double target, double lim, int param, int type) {

  const int n = lprob.size();
  NumericVector lprobFl(n);
  NumericVector prob;
  bool out = false;

  for (int i = 0; i < n; ++i) {
    lprobFl[i] = lprob[i] + lim;
  }

  if (param == 0) {

    prob = Rcpp::plogis(lprobFl);

  } else {

    prob = Rcpp::plogis(lprobFl) - Rcpp::plogis(lprob);

  }

  if (type == 0) {

    out = mean(prob) < target;

  } else  {

    out = mean(prob) > target;

  }

  return out;

}


// [[Rcpp::export(find_b0_cpp)]]
double find_b0(NumericVector lprob, double target, double tol, int maxIter = 100, double intLen = 0.001) {

  double ll = -10;
  double ul = 10;

  while (!check_limit(lprob, target, ll, 0, 0)) {

    ll = ll * 2;

  }

  while (!check_limit(lprob, target, ul, 0, 1)) {

    ul = ul * 2;

  }


  const int n = lprob.size();
  NumericVector lprobFull(n);
  double b0, prev = target + 2 * tol;
  int counter = 0;

  while (std::abs(prev - target) > tol) {

    if (counter == maxIter) {
      Rcpp::stop("Maximum number of iterations reached. Adjust interval limits.");
    }

    b0 = (ll + ul) * 0.5;

    for (int i = 0; i < n; ++i) {
      lprobFull[i] = lprob[i] + b0;
    }

    NumericVector prob = Rcpp::plogis(lprobFull);
    prev = mean(prob);

    if (prev < target) {
      ll = b0;
    } else {
      ul = b0;
    }

    ++counter;
  }

  return b0;

}

// [[Rcpp::export(find_bz_cpp)]]
double find_bz(NumericVector lprob, double target, double tol, int maxIter = 100, double intLen = 0.001) {

  double ll = -10;
  double ul = 10;

  while (!check_limit(lprob, target, ll, 1, 0)) {

    ll = ll * 2;

  }

  while (!check_limit(lprob, target, ul, 1, 1)) {

    ul = ul * 2;

  }

  const int n = lprob.size();
  NumericVector lprobFull(n);
  double bz, rd = target + 2 * tol;
  int counter = 0;

  while (std::abs(rd - target) > tol) {

    if (counter == maxIter) {
      Rcpp::stop("Maximum number of iterations reached. Adjust interval limits.");
    }

    bz = (ll + ul) * 0.5;

    for (int i = 0; i < n; ++i) {
      lprobFull[i] = lprob[i] + bz;
    }

    NumericVector diff = Rcpp::plogis(lprobFull) - Rcpp::plogis(lprob);
    rd = mean(diff);

    if (rd < target) {
      ll = bz;
    } else {
      ul = bz;
    }

    ++counter;
  }

  return bz;

}

