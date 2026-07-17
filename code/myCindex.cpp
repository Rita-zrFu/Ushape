// calculates my C-index, original and smooth versions
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double myCindex(NumericVector H, NumericVector time, NumericVector delta) {
  int n = H.size();
  double concordant = 0.0, tied = 0.0, comparable = 0.0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      if (time[j] > time[i] && delta[i] == 1) {
        comparable += 1.0;
        if (H[j] < H[i]) concordant += 1.0;
        else if (H[i] == H[j]) tied += 1.0;
      }
    }
  }
  return (concordant + tied / 2.0) / comparable;
}

// [[Rcpp::export]]
double myCindex_smooth(NumericVector H, NumericVector time,
                       NumericVector delta, double sigma_h) {
  // smoothed version: replace indicator 1{H_i > H_j} by sigmoid( - (H_i - H_j) * sigma_h )
  int n = H.size();
  double num1 = 0.0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (time[i] > time[j] && delta[j] == 1) {
        double sigmoid = 1.0 / (1.0 + std::exp(sigma_h * (H[i] - H[j])));
        num1 += sigmoid;
      }
    }
  }
  return 2.0 * num1 / (n * (n - 1.0));
}
