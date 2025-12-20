#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm> 
#include <random>

using namespace Rcpp;

double inner_compute(double cross_inner, 
                     double self_inner1, 
                     double self_inner2) {
  double value = 0.0;
  double denominator = (1.0 + self_inner1) * (1.0 + self_inner2);
  denominator = std::sqrt(denominator);
  value = (1.0 + cross_inner) / denominator; 
  value = std::asin(value);
  return value;
}

double compute_ph_value(Eigen::MatrixXd &asin_X,
                        Eigen::MatrixXd &asin_Y, 
                        Eigen::VectorXd &asin_X_col, 
                        Eigen::VectorXd &asin_Y_col, 
                        double asin_X_mean, 
                        double asin_Y_mean, 
                        int num) {
  double tmp_A, tmp_B;
  double ph_value = 0.0; 
  for (int i = 0; i < num; i++) {
    for (int j = (i + 1); j < num; j++) {
      tmp_A = asin_X(i, j) - asin_X_col(i) - asin_X_col(j) + asin_X_mean;
      tmp_B = asin_Y(i, j) - asin_Y_col(i) - asin_Y_col(j) + asin_Y_mean;
      ph_value += 2.0 * tmp_A * tmp_B;
    }
  }
  
  for (int i = 0; i < num; i++) {
    tmp_A = asin_X(i, i) - 2 * asin_X_col(i) + asin_X_mean;
    tmp_B = asin_Y(i, i) - 2 * asin_Y_col(i) + asin_Y_mean;
    ph_value += tmp_A  * tmp_B;
  }
  return ph_value;
}

Eigen::MatrixXd shuffle_matrix(Eigen::MatrixXd A, 
                               std::vector<int> index, 
                               int num) {
  Eigen::MatrixXd A_tmp = Eigen::MatrixXd::Zero(num, num);
  for (int i = 0; i < num; i++) {
    for (int j = (i + 1); j < num; j++) {
      A_tmp(i, j) = A(index[i], index[j]);
      A_tmp(j, i) = A_tmp(i, j); 
    }
  }
  
  for (int i = 0; i < num; i++) {
    A_tmp(i, i) = A(index[i], index[i]);
  }
  return A_tmp;
}

Eigen::MatrixXd shuffle_vector(Eigen::VectorXd vec, 
                               std::vector<int> index, 
                               int num) {
  Eigen::VectorXd vec_tmp = Eigen::VectorXd::Zero(num);
  
  for (int i = 0; i < num; i++) {
    vec_tmp(i) = vec(index[i]);
  }
  return vec_tmp;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List rcpp_pht(Eigen::MatrixXd &X, Eigen::MatrixXd &Y, int R) {
  double ph_value = 0.0; 
  int num = X.rows();
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(num, num);
  Eigen::MatrixXd B = Eigen::MatrixXd::Identity(num, num);
  A = A * (PI / 2);
  B = B * (PI / 2);
  
  
  Eigen::MatrixXd cov_X = X * X.transpose();
  Eigen::MatrixXd cov_Y = Y * Y.transpose();
  
  for (int i = 0; i < num; i++) {
    for (int j = (i + 1); j < num; j++) {
      A(i, j) = inner_compute(cov_X(i, j), cov_X(i, i), cov_X(j, j));
      A(j, i) = A(i, j); 
      B(i, j) = inner_compute(cov_Y(i, j), cov_Y(i, i), cov_Y(j, j));
      B(j, i) = B(i, j); 
    }
  }
  // std::cout << "Matrix A: " << std::endl;
  // std::cout << A << std::endl;
  // std::cout << "Matrix B: " << std::endl;
  // std::cout << B << std::endl;
  
  Eigen::VectorXd A_col_mean = A.colwise().mean();
  Eigen::VectorXd B_col_mean = B.colwise().mean();
  double A_mean = A_col_mean.mean();
  double B_mean = B_col_mean.mean();
  
  
  // compute test statistic:
  ph_value = compute_ph_value(A, B, 
                              A_col_mean, B_col_mean, 
                              A_mean, B_mean, 
                              num);
  
  // compute permuted statistic:
  double p_value = 1.0;
  std::vector<int> random_index(num);
  for (int i = 0; i < num; i++) {
    random_index[i] = i;
  }
  
  Eigen::VectorXd permuted_ph_value = Eigen::VectorXd::Zero(R);
  std::random_device rd;
  std::mt19937 g(rd());
  for (int r = 0; r < R; r++) {
    std::shuffle(random_index.begin(), random_index.end(), g);
    A = shuffle_matrix(A, random_index, num);
    A_col_mean = shuffle_vector(A_col_mean, random_index, num);
    permuted_ph_value(r) = compute_ph_value(A, B, 
                                            A_col_mean, 
                                            B_col_mean, 
                                            A_mean, B_mean, 
                                            num);
    if (permuted_ph_value(r) > ph_value) {
      p_value += 1.0;
    }
  }
  p_value /= (1.0 + R);
  
  List res = List::create(Named("statistic") = ph_value, 
                     Named("permuted.statistic") = permuted_ph_value, 
                     Named("p.value") = p_value); 
    
  return res; 
}

