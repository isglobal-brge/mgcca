#ifndef _rfunctions_UTILS_H
#define _rfunctions_UTILS_H


//
// REFERENCES:
//  Function from:
//    https://github.com/jaredhuling/rfunctions
//


#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cmath>


// using Eigen::MatrixXd;
// using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Lower;



//computes X'WX where W is diagonal (input w as vector)
Eigen::MatrixXd xtwx(const Eigen::MatrixXd& xx, const Eigen::MatrixXd& ww);

//computes X'SX where S is not diagonal (input ss as matrix)
Eigen::MatrixXd xtsx(const Eigen::MatrixXd& xx, const Eigen::MatrixXd& ss);

//computes X'X
Eigen::MatrixXd xtx(const Eigen::MatrixXd& xx);

//computes XX'
Eigen::MatrixXd xxt(const Eigen::MatrixXd& xx);

//solve Ax = b using conjugate gradient
Eigen::MatrixXd conjugate_gradient(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, int maxit, double tol);

#endif
