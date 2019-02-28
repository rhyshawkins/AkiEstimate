#ifndef generalisedeigenproblem_hpp
#define generalisedeigenproblem_hpp

#include "spec1dmatrix.hpp"

//
// Copied from
//
//   http://eigen.tuxfamily.org/index.php?title=Lapack
//

extern "C" void dggev_(const char* JOBVL, const char* JOBVR, const int* N,
                       const double* A, const int* LDA, const double* B, const int* LDB,
                       double* ALPHAR, double* ALPHAI, double* BETA,
                       double* VL, const int* LDVL, double* VR, const int* LDVR,
                       double* WORK, const int* LWORK, int* INFO);

// Generalised Eigen-Problem
// Solve:
// A * v(j) = lambda(j) * B * v(j).
//
// v are the eigenvectors and are stored in v.
// lambda are the eigenvalues and are stored in lambda.
// The eigenvalues are stored as: (lambda(:, 1) + lambda(:, 2)*i)./lambda(:, 3)
//
// returns true on success.
template
<
  typename real
>
bool GEP(const Spec1DMatrix<real> &_A,
	 const Spec1DMatrix<real> &_B,
	 Spec1DMatrix<real> &_work,
	 Spec1DMatrix<real> &_u,
	 Spec1DMatrix<real> &_v,
	 Spec1DMatrix<real> &_lambda)
{
  Spec1DMatrix<double> A;
  Spec1DMatrix<double> B;
  Spec1DMatrix<double> work;
  Spec1DMatrix<double> u;
  Spec1DMatrix<double> v;
  Spec1DMatrix<double> lambda;

  int N = _A.cols();
  if (_B.cols() != N || _A.rows()!=N || _B.rows() != N) {
    return false;
  }

  A.resize(N,N);
  B.resize(N,N);
  v.resize(N,N);
  lambda.resize(N,3);

  for (int j = 0; j < N; j ++) {
    for (int i = 0; i < N; i ++) {
      A(i, j) = value(_A(i, j));
      B(i, j) = value(_B(i, j));
    }
  }

  if (!GEP<double>(A, B, work, u, v, lambda)) {
    return false;
  }

  _u.resize(N, N);
  _v.resize(N, N);
  _lambda.resize(N, 3);

  for (int j = 0; j < N; j ++) {
    for (int i = 0; i < N; i ++) {
      _u(i, j) = u(i, j);
      _v(i, j) = v(i, j);
    }

    _lambda(j, 0) = lambda(j, 0);
    _lambda(j, 1) = lambda(j, 1);
    _lambda(j, 2) = lambda(j, 2);
  }
  
  return true;
}

template
<>
bool GEP<double>(const Spec1DMatrix<double> &A,
		 const Spec1DMatrix<double> &B,
		 Spec1DMatrix<double> &work,
		 Spec1DMatrix<double> &u,
		 Spec1DMatrix<double> &v,
		 Spec1DMatrix<double> &lambda)
{
  int N = A.cols(); 
  if (B.cols() != N || A.rows()!=N || B.rows() != N) {
    return false;
  }

  u.resize(N, N);
  v.resize(N, N);
  lambda.resize(N, 3);

  int LDA = N;
  int LDB = N;
  int LDV = N;

  double WORKDUMMY;
  int LWORK = -1; // Request optimum work size.
  int INFO = 0;
  
  double *alphar = lambda.col(0);
  double *alphai = lambda.col(1);
  double *beta   = lambda.col(2);

  // Get the optimum work size.
  dggev_("V", "V", &N, A.data(), &LDA, B.data(), &LDB, alphar, alphai, beta,
	 u.data(), &LDV, v.data(), &LDV, &WORKDUMMY, &LWORK, &INFO);

  LWORK = int(WORKDUMMY) + 32;
  work.resize(LWORK, 1);

  dggev_("V", "V", &N, A.data(), &LDA, B.data(), &LDB, alphar, alphai, beta,
	 u.data(), &LDV, v.data(), &LDV, work.data(), &LWORK, &INFO);

  return INFO==0;
}

#endif // generalisedeigenproblem_hpp
