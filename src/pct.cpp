//
//  pct.cpp
//  pct
//
//  Created by wangtingyin on 2021/2/26.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>
#include <cerrno>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

double inner_compute_pc(double cross_inner,
                        double self_inner1,
                        double self_inner2){
    double value = 0.0;
    double denominator = self_inner1 * self_inner2;
    denominator = std::sqrt(denominator);
    value = cross_inner / denominator;
    
    if (denominator == 0.0) {
      value = 0.0;
    }
    value = std::fabs(value) >= 1.0 ? 0.0 : std::acos(value);
//std::cout << value << endl;
    // errno = 0;
    // value = std::acos(value);
    // if (errno == EDOM) {
    //   value = 0.0;
    // }
    return value;
}

double compute_pc_value(std::vector<std::vector<std::vector<double> > > &asin_X,
                        std::vector<std::vector<std::vector<double> > > &asin_Y,
                        Eigen::MatrixXd &asin_X_col,
                        Eigen::MatrixXd &asin_Y_col,
                        Eigen::VectorXd &asin_X_mat,
                        Eigen::VectorXd &asin_Y_mat,
                        int num, bool term_S2){
    double tmp_A, tmp_B;
    double pc_value = 0.0;
    double stat;
    for(int r = 0; r < num; r++) {
        for(int k = 0; k < num; k++) {
            for(int l = (k+1); l < num; l++) {
                tmp_A = asin_X[k][l][r] - asin_X_col(k, r) - asin_X_col(l, r) + asin_X_mat(r);
                tmp_B = asin_Y[k][l][r] - asin_Y_col(k, r) - asin_Y_col(l, r) + asin_Y_mat(r);
                pc_value += 2.0 * tmp_A * tmp_B;
            }
        }
    }
    
    for(int r = 0; r < num; r++){
        for(int k = 0; k < num; k++){
            tmp_A = asin_X[k][k][r] - 2*asin_X_col(k,r) + asin_X_mat(r);
            tmp_B = asin_Y[k][k][r] - 2*asin_Y_col(k,r) + asin_Y_mat(r);
            pc_value += tmp_A * tmp_B;
        }
    }
    
    pc_value = pc_value / (num*num*num);
    stat = pc_value;
    
    if (term_S2) {
      double S_2 = 0.0;
      for(int i = 0; i < num; i++){
        for(int j = 0; j < num; j++){
          for(int k = 0; k < num; k++){
            for(int l = 0; l < num; l++){
              for(int r = 0; r < num; r++){
                S_2 += asin_X[i][l][k] * asin_Y[j][r][k];
              }
            }
          }
        }
      }
      S_2 = S_2 / (num * num * num * num * num);
      stat = num * pc_value *pc_value / (9.869604 - S_2);
    }
    
    return stat;
}

std::vector<std::vector<std::vector<double> > >  shuffle_3D_array_pc(std::vector<std::vector<std::vector<double> > > &A,
                                                                     std::vector<int> index,
                                                                     int num){
    std::vector<std::vector<std::vector<double> > > A_tmp (num, std::vector<std::vector<double> >(num, std::vector <double>(num, 0.0)));
    for(int r = 0; r < num; r++){
        for(int k = 0; k < num; k++){
            for(int l = (k+1); l < num; l++){
                A_tmp[k][l][r] = A[index[k]][index[l]][index[r]];
                A_tmp[l][k][r] = A_tmp[k][l][r];
            }
        }
    }
    
    for (int r = 0; r < num; r++){
        for(int k = 0; k < num; k++){
            A_tmp[k][k][r] = A[index[k]][index[k]][index[r]];
        }
    }
    return A_tmp;
}

Eigen::MatrixXd shuffle_matrix_pc(Eigen::MatrixXd A,
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

Eigen::VectorXd shuffle_vector_pc(Eigen::VectorXd vec,
                                  std::vector<int> index,
                                  int num) {
  Eigen::VectorXd vec_tmp = Eigen::VectorXd::Zero(num);
  
  for (int i = 0; i < num; i++) {
    vec_tmp(i) = vec(index[i]);
  }
  return vec_tmp;
}

double scalar_product_pc(Eigen::MatrixXd X, int k, int l, int r){
    double product = 0.0;
    // int num = X.rows();
    // for (int i = 0; i < num; i++){
    //    product = product + (X(i,k)-X(i,r))*(X(i,l)-X(i,r));
    // }
    product = (X.row(k) - X.row(r)).dot(X.row(l) - X.row(r));
    return product;
}

double mean_1_dim_pc(std::vector<std::vector<std::vector<double> > > &A, int num, int l, int r){
    double sum = 0.0;
    double mean;
    for(int k = 0; k < num; k++){
        sum = sum + A[k][l][r];
    }
    mean = sum / num;
    return mean;
}

double mean_2_dim_pc(std::vector<std::vector<std::vector<double> > > &A, int num, int r){
    double sum = 0.0;
    double mean;
    for(int k = 0; k < num; k++) {
        for(int l = 0; l < num; l++) {
            sum = sum + A[k][l][r];
        }
    }
    mean = sum / (num * num);
    return mean;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List rcpp_pcov(Eigen::MatrixXd &X, Eigen::MatrixXd &Y, int R, bool term_S2){
    double pcov_value = 0.0;
    int num = X.rows();
    std::vector<std::vector<std::vector<double> > > A (num, std::vector<std::vector<double> >(num, std::vector<double>(num, 0.0)));
    std::vector<std::vector<std::vector<double> > > B (num, std::vector<std::vector<double> >(num, std::vector<double>(num, 0.0)));
    Eigen::MatrixXd mean_mat_A = Eigen::MatrixXd::Zero(num, num);
    Eigen::MatrixXd mean_mat_B = Eigen::MatrixXd::Zero(num, num);
    Eigen::VectorXd mean_col_A = Eigen::VectorXd::Zero(num);
    Eigen::VectorXd mean_col_B = Eigen::VectorXd::Zero(num);
    
    for(int r = 0; r < num; r++){
        for(int k = 0; k < num; k++){
            for(int l = (k + 1); l < num; l++) {
                A[k][l][r] = inner_compute_pc(scalar_product_pc(X, k, l, r),
                                              scalar_product_pc(X, k, k, r),
                                              scalar_product_pc(X, l, l, r));
                A[l][k][r] = A[k][l][r];
                B[k][l][r] = inner_compute_pc(scalar_product_pc(Y, k, l, r), 
                                              scalar_product_pc(Y, k, k, r), 
                                              scalar_product_pc(Y, l, l, r));
                B[l][k][r] = B[k][l][r];
            }
            A[k][k][r] = 0.0;
            B[k][k][r] = 0.0;
        }
        for(int k = 0; k < num; k++){
          A[k][r][r] = 0.0;
          B[k][r][r] = 0.0;
        }
        for(int k = 0; k < num; k++){
          A[r][k][r] = 0.0;
          B[r][k][r] = 0.0;
        }
    }
    
    // for(int r = 0; r < num; r++) {
    //   std::cout << "=================" << std::endl;
    //   for(int k = 0; k < num; k++) {
    //     for(int l = 0; l < num; l++) { 
    //       std::cout << A[k][l][r] << " ";
    //     }
    //     std::cout << std::endl;
    //   }
    // }
    
    for(int l = 0; l < num; l++){
        for(int r = 0; r < num; r++){
            mean_mat_A(l, r) = mean_1_dim_pc(A, num, l, r);
            mean_mat_B(l, r) = mean_1_dim_pc(B, num, l, r);
        }
    }
    
    for(int r = 0; r < num; r++){
        mean_col_A[r] = mean_2_dim_pc(A, num, r);
        mean_col_B[r] = mean_2_dim_pc(B, num, r);
    }
    // compute test statistic:
    // std::cout << "compute statistic" << std::endl;
    pcov_value =  compute_pc_value(A,
                                   B,
                                   mean_mat_A,
                                   mean_mat_B,
                                   mean_col_A,
                                   mean_col_B,
                                   num, term_S2);
    
    // compute permuted statistic:
    double p_value = 1.0;
    std::vector<int> random_index(num);
    for (int i = 0; i < num; i++) {
      random_index[i] = i;
    }
    
    Eigen::VectorXd permuted_pcov_value = Eigen::VectorXd::Zero(R);
    std::random_device rd;
    std::mt19937 g(rd());
    for (int r = 0; r < R; r++) {
      std::shuffle(random_index.begin(), random_index.end(), g);
      std::vector<std::vector<std::vector<double> > > As;
      Eigen::MatrixXd mean_mat_As;
      Eigen::VectorXd mean_col_As;
      As = shuffle_3D_array_pc(A, random_index, num);
      mean_mat_As = shuffle_matrix_pc(mean_mat_A, random_index, num);
      mean_col_As = shuffle_vector_pc(mean_col_A, random_index, num);
      permuted_pcov_value(r) = compute_pc_value(As, B,
                                                mean_mat_As,
                                                mean_mat_B,
                                                mean_col_As,
                                                mean_col_B,
                                                num, term_S2);
      if (permuted_pcov_value(r) > pcov_value) {
        p_value += 1.0;
      }
    }
    p_value /= (1.0 + R);
    
    List res = List::create(Named("statistic") = pcov_value,
                            Named("permuted.statistic") = permuted_pcov_value,
                            Named("p.value") = p_value);
      
    return res;
    
}


