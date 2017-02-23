transCpp <- 'using Eigen::Map;
using Eigen::MatrixXi;
// Map the integer matrix AA from R
const Map <MatrixXi> A(as<Map<MatrixXi>  >(AA));
// evaluate and return the transpose of A
const MatrixXi      At(A.transpose());
return wrap(At);'

ftrans <- cxxfunction(signature(AA = "matrix"), transCpp,
                      plugin = "RcppEigen")

prodCpp <- 'typedef Eigen::Map<Eigen::MatrixXi>  MapMati;
const MapMati    B(as<MapMati>(BB));
const MapMati    C(as<MapMati>(CC));
return List::create(Named("B %*% C")         = B * C,
                    Named("crossprod(B, C)") = B.adjoint() * C);'

fprod <- cxxfunction(signature(BB = "matrix", CC = "matrix"),
                     prodCpp, "RcppEigen")

crossprodCpp <- 'using Eigen::Map;
using Eigen::MatrixXi;
using Eigen::Lower;
const Map <MatrixXi> A(as<Map<MatrixXi>  >(AA));
const int           m(A.rows()), n(A.cols());
MatrixXi            AtA(MatrixXi(n, n).setZero().
                        selfadjointView<Lower>().rankUpdate(A.adjoint()));
MatrixXi            AAt(MatrixXi(m, m).setZero().
                        selfadjointView<Lower>().rankUpdate(A));
return List::create(Named("crossprod(A)")  = AtA, 
                    Named("tcrossprod(A)") = AAt);'

fcprd <- cxxfunction(signature(AA = "matrix"), crossprodCpp, "RcppEigen")



