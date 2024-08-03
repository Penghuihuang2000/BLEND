#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List compute_U_n(NumericVector X_n, List phi, List psi_n_mu_n, List U_n, NumericVector cumsum_Ms, NumericVector cumsum_Ms_Ms, int num_ct) {
  int G = X_n.size(); // number of genes
  int CT = phi.size(); // number of cell types
  List U_n_new = clone(U_n);
  List phi_g(CT);

  for (int g = 0; g < G; ++g) {

    for (int ct = 0; ct < CT; ++ct) {
      phi_g[ct] = as<NumericMatrix>(phi[ct])(g, _);
    }
    List prob_list = clone(psi_n_mu_n);
    for (int ct = 0; ct < CT; ++ct) {
      prob_list[ct] = as<NumericVector>(prob_list[ct]) * as<NumericVector>(phi_g[ct]);
    }

    double prob_sum = 0.0;
    for (int ct = 0; ct < CT; ++ct) {
      NumericVector prob_vec = prob_list[ct];
      for (int i = 0; i < prob_vec.size(); ++i) {
        prob_sum += prob_vec[i];
      }
    }

    NumericVector U_n_sample(num_ct);
    int index = 0;
    for (int ct = 0; ct < CT; ++ct) {
      NumericVector prob_vec = prob_list[ct];
      for (int i = 0; i < prob_vec.size(); ++i) {
        U_n_sample[index] = X_n[g] * prob_vec[i] / prob_sum;
        ++index;
      }
    }

    for (int ct = 0; ct < CT; ++ct) {
      int start_idx = cumsum_Ms_Ms[ct] - 1;
      int end_idx = cumsum_Ms[ct] - 1;
      as<NumericMatrix>(U_n_new[ct])(g, _) = U_n_sample[Range(start_idx, end_idx)];
    }
  }

  return U_n_new;
}
