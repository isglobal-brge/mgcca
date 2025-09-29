

#include "utils.h"


//computes X'WX where W is diagonal (input w as vector)
Eigen::MatrixXd xtwx(const Eigen::MatrixXd& xx, const Eigen::MatrixXd& ww) {
	const int n(xx.cols());
	Eigen::MatrixXd AtWA(Eigen::MatrixXd(n, n).setZero().
			 selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.asDiagonal()));
	return (AtWA);
}

//computes X'SX where S is not diagonal (input ss as matrix)
Eigen::MatrixXd xtsx(const Eigen::MatrixXd& xx, const Eigen::MatrixXd& ss) {
	const int n(xx.cols());
	Eigen::MatrixXd AtSA(Eigen::MatrixXd(n, n).setZero().
							 selfadjointView<Lower>().rankUpdate(xx.adjoint() * ss));
	return (AtSA);
}

Eigen::MatrixXd xtx(const Eigen::MatrixXd& xx) {
	const int n(xx.cols());
	Eigen::MatrixXd AtA(Eigen::MatrixXd(n, n).setZero().
							selfadjointView<Lower>().rankUpdate(xx.adjoint()));
	return (AtA);
}

Eigen::MatrixXd xxt(const Eigen::MatrixXd& xx) {
	const int m(xx.rows());
	Eigen::MatrixXd AtA(Eigen::MatrixXd(m, m).setZero().
		selfadjointView<Lower>().rankUpdate(xx));
	return (AtA);
}


Eigen::MatrixXd conjugate_gradient(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, int maxit, double tol)
{

	const int n(A.cols());
	Eigen::VectorXd x(n);
	Eigen::VectorXd r(n);
	Eigen::VectorXd p(n);
	Eigen::VectorXd Ap(n);
	x.fill(0);

	double rsold;
	double rsnew;
	double alpha;
	int iters = maxit;

	r = b;
	p = r;
	rsold = r.squaredNorm();

	for (int i = 0; i < maxit; i++) {
		Ap = A * p;
		alpha = rsold / (p.transpose() * Ap);
		x = x + (alpha * p.array()).matrix();
		r = r - (alpha * Ap.array()).matrix();
		rsnew = r.squaredNorm();
		if (sqrt(rsnew) < tol) {
			iters = i;
			break;
		}
		p = r + ((rsnew / rsold) * p.array()).matrix();
		rsold = rsnew;
	}
	return(x);
}


