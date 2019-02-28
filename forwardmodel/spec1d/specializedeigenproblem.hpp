#pragma once
#ifndef specializedeigenproblem_hpp
#define specializedeigenproblem_hpp

#include <math.h>

#include "spec1dmatrix.hpp"

extern "C" {
  void dgghrd_(char *compq,
	       char *compz,
	       int *n,
	       int *ilo,
	       int *ihi,
	       double *a,
	       int *lda,
	       double *b,
	       int *ldb,
	       double *q,
	       int *ldq,
	       double *z,
	       int *ldz,
	       int *info);

  void dhgeqz_(char *job,
	       char *compq,
	       char *compz,
	       int *n,
	       int *ilo,
	       int *ihi,
	       double *h,
	       int *ldh,
	       double *t,
	       int *ldt,
	       double *alphar,
	       double *alphai,
	       double *beta,
	       double *q,
	       int *ldq,
	       double *z,
	       int *ldz,
	       double *work,
	       int *lwork,
	       int *info);

  void dtgevc_(char *side,
	       char *howmny,
	       int *select,
	       int *n,
	       double *s,
	       int *lds,
	       double *p,
	       int *ldp,
	       double *vl,
	       int *ldvl,
	       double *vr,
	       int *ldvr,
	       int *mm,
	       int *m,
	       double *work,
	       int *info);
}

template
<
  typename real
>
bool SpecializedEigenProblem(Spec1DMatrix<real> &_A,
			     Spec1DMatrix<real> &_B,
			     Spec1DMatrix<real> &_work,
			     Spec1DMatrix<real> &_u,
			     Spec1DMatrix<real> &_v,
			     Spec1DMatrix<real> &_lambda,
			     Spec1DMatrix<real> &_Q,
			     Spec1DMatrix<real> &_Z)
{
  fprintf(stderr, "SpecializedEigenProblem: unimplemented\n");
  return false;
}

template
<>
bool SpecializedEigenProblem<double>(Spec1DMatrix<double> &A,
				     Spec1DMatrix<double> &B,
				     Spec1DMatrix<double> &work,
				     Spec1DMatrix<double> &VL,
				     Spec1DMatrix<double> &VR,
				     Spec1DMatrix<double> &lambda,
				     Spec1DMatrix<double> &Q,
				     Spec1DMatrix<double> &Z)
{
  int N = A.rows();
  int lwork = 6*N; // DTGEVC requires workspace 6*N which is largest

  work.resize(6*N, 1);
  lambda.resize(N, 3);

  //
  // Initialize eigen vectors to Identity (this resizes as well)
  //
  VL.setIdentity(N);
  VR.setIdentity(N);

  //
  // DGGHRD : Reduce A,B to General Hessenberg form.
  //
  int ilo = 1;
  int ihi = N;

  int lda = N;
  int ldb = N;
  int ldq = N;
  int ldz = N;

  char compI[] = "I";
  int info;
  
  dgghrd_(compI,
	  compI,
	  &N,
	  &ilo,
	  &ihi,
	  A.data(),
	  &lda,
	  B.data(),
	  &ldb,
	  VL.data(),
	  &ldq,
	  VR.data(),
	  &ldz,
	  &info);

  if (info != 0) {
    ERROR("Failed to reduce to general Hessenberg form");
    return false;
  }

  //
  // DHGEQZ : Compute Eigen values using Schur decomposition
  //
  char schur[] = "S";
  char compV[] = "V";
  
  dhgeqz_(schur,
	  compV,
	  compV,
	  &N,
	  &ilo,
	  &ihi,
	  A.data(),
	  &lda,
	  B.data(),
	  &ldb,
	  lambda.col(0),
	  lambda.col(1),
	  lambda.col(2),
	  VL.data(),
	  &ldq,
	  VR.data(),
	  &ldz,
	  work.data(),
	  &lwork,
	  &info);
 
  //
  // Copy VL, VR to Q, Z (assignment operator resizes)
  //
  Q = VL;
  Z = VR;
  
  //
  // DTGEVC : Compute Eigen vectors
  //
  char both[] = "B";
  char allback[] = "B";
  int mm = N;
  int m;
  
  dtgevc_(both,
	  allback,
	  nullptr,
	  &N,
	  A.data(),
	  &lda,
	  B.data(),
	  &ldb,
	  VL.data(),
	  &ldq,
	  VR.data(),
	  &ldz,
	  &mm,
	  &m,
	  work.data(),
	  &info);

  if (info != 0) {
    ERROR("Failed to compute eigen vectors\n");
    return false;
  }

  return true;
}


//
// Based upon 2x2 real solver of DLALN2 Lapack subroutine
//
bool
Solve2x2(double a11,
	 double a12,
	 double a21,
	 double a22,
	 double b1,
	 double b2,
	 double &x1,
	 double &x2)
{
  double ur11;
  double cr21;
  double ur12;
  double cr22;

  bool rswap;
  bool xswap;

  //
  // Pivot to place largest value at 11
  //
  if (fabs(a11) > fabs(a12)) {
    if (fabs(a11) > fabs(a21)) {
      if (fabs(a11) > fabs(a22)) {
	// 0
	ur11 = a11;
	cr21 = a21;
	ur12 = a12;
	cr22 = a22;

	rswap = false;
	xswap = false;
	
      } else {
	// 3
	ur11 = a22;
	cr21 = a12;
	ur12 = a21;
	cr22 = a11;

	rswap = true;
	xswap = true;
	
      }
    } else {
      if (fabs(a21) > fabs(a22)) {
	// 2
	ur11 = a21;
	cr21 = a11;
	ur12 = a22;
	cr22 = a12;

	rswap = true;
	xswap = false;
	
      } else {
	// 3
	ur11 = a22;
	cr21 = a12;
	ur12 = a21;
	cr22 = a11;

	rswap = true;
	xswap = true;
      }
    }
  } else {
    if (fabs(a12) > fabs(a21)) {
      if (fabs(a12) > fabs(a22)) {
	// 1
	ur11 = a12;
	cr21 = a22;
	ur12 = a11;
	cr22 = a21;

	rswap = false;
	xswap = true;
	
      } else {
	// 3
	ur11 = a22;
	cr21 = a12;
	ur12 = a21;
	cr22 = a11;

	rswap = true;
	xswap = true;
      }
    } else {
      if (fabs(a21) > fabs(a22)) {
	// 2 
	ur11 = a21;
	cr21 = a11;
	ur12 = a22;
	cr22 = a12;

	rswap = true;
	xswap = false;
      } else {
	// 3
	ur11 = a22;
	cr21 = a12;
	ur12 = a21;
	cr22 = a11;

	rswap = true;
	xswap = true;
      }
    }
  }

  double ur11r = 1.0/ur11;
  double lr21 = ur11r*cr21;
  double ur22 = cr22 - ur12*lr21;

  double br1, br2;
  if (rswap) {
    br1 = b2;
    br2 = b1;
  } else {
    br1 = b1;
    br2 = b2;
  }

  br2 -= lr21*br1;

  double xr2 = br2/ur22;
  double xr1 = br1*ur11r - xr2*(ur11r*ur12);
  
  if (xswap) {
    x1 = xr2;
    x2 = xr1;
  } else {
    x1 = xr1;
    x2 = xr2;
  }
    
  return true;
}
template
<
  typename real
>
bool SpecializedEigenProblemAdjoint(Spec1DMatrix<real> &S,
				    Spec1DMatrix<real> &P,
				    Spec1DMatrix<real> &Q,
				    Spec1DMatrix<real> &Z,
				    real alpha,
				    Spec1DMatrix<real> &B,
				    Spec1DMatrix<real> &lambda,
				    Spec1DMatrix<real> &work)
{
  ERROR("Unimplemented");
  return false;
}

void sepdump(const char *filename, Spec1DMatrix<double> &m)
{
  FILE *fp = fopen(filename, "w");
  int N = m.rows();
  int M = m.cols();
  for (int i = 0; i < N; i ++) {
    for (int j = 0; j < M; j ++) {
      
      fprintf(fp, "%25.19e ", m(i, j));
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
	     
	  
template
<>
bool SpecializedEigenProblemAdjoint<double>(Spec1DMatrix<double> &S,
					    Spec1DMatrix<double> &P,
					    Spec1DMatrix<double> &Q,
					    Spec1DMatrix<double> &Z,
					    double alpha,
					    Spec1DMatrix<double> &B,
					    Spec1DMatrix<double> &lambda,
					    Spec1DMatrix<double> &work)
{
  int N = S.rows();
  work.resize(2*N, 1);

  if (B.rows() != N || lambda.rows() != N) {
    FATAL("B/Lambda no. rows incorrect");
    return false;
  }

  if (B.cols() != lambda.cols()) {
    FATAL("B/Lambda columns mismatch\n");
    return false;
  }
  

  for (int k = 0; k < B.cols(); k ++) {
  
    //
    // Multiply B by Z^T and store in work[0 .. N - 1]
    //
    for (int i = 0; i < N; i ++) {
      double s = 0.0;
      for (int j = 0; j < N; j ++) {
	s += Z(j, i) * B(j, k);
      }
      work(i, 0) = s;
    }

    //
    // Forward substitute S^T - alpha P^T to solve for Q^T lambda and store in work[N .. 2*N - 1]
    //
    for (int i = 0; i < N; i ++) {

      if (i < (N - 1) && S(i + 1, i) != 0.0) {

	//
	// 2x2
	//
	double a = 0.0;
	double c = 0.0;
	for (int j = 0; j < i; j ++) {
	  a += (S(j, i) - alpha*P(j, i)) * work(N + j, 0);
	  c += (S(j, i + 1) - alpha*P(j, i + 1)) * work(N + j, 0);
	}

	double bi = S(i, i) - alpha*P(i, i);
	double bi1 = S(i + 1, i) - alpha*P(i + 1, i);

	double di = S(i, i + 1) - alpha*P(i, i + 1);
	double di1 = S(i + 1, i + 1) - alpha*P(i + 1, i + 1);
	
	
	// a + bi Q^T lambda_i + bi1 Q^T lambda_{i + 1} = work[i]
	// c + di Q^T lambda_i + di1 Q^T lambda_{i + 1} = work[i + 1]
	Solve2x2(bi, bi1, di, di1, work(i, 0) - a, work(i + 1, 0) - c, work(N + i, 0), work(N + i + 1, 0));
	i ++;
	
      } else {

	//
	// 1x1
	//
	double denom = (S(i, i) - alpha*P(i, i));
	double s = 0.0;
	if (denom == 0.0) {
	  work(N + i, 0) = 0.0;
	} else {
	  
	  for (int j = 0; j < i; j ++) {
	    s += (S(j, i) - alpha*P(j, i)) * work(N + j, 0);
	  }
	  
	  work(N + i, 0) = (work(i, 0) - s)/denom;
	}

      }
    }
    
    //
    // Multiply work[N .. 2*N - 1] by Q and store in lambda
    //
    for (int i = 0; i < N; i ++) {
      double s = 0.0;
      for (int j = 0; j < N; j ++) {
	s += Q(i, j) * work(N + j, 0);
      }
      lambda(i, k) = s;
    }
  }

  return true;
}

void
SpecializedEigenProblemVerify(Spec1DMatrix<double> &S,
			      Spec1DMatrix<double> &Q,
			      Spec1DMatrix<double> &Z,
			      Spec1DMatrix<double> &A,
			      Spec1DMatrix<double> &work)
{
  int N = S.rows();
  
  work.resize(N, N);
  A.resize(N, N);

  //
  // Compute Q S
  //
  for (int i = 0; i < N; i ++) {
    
    for (int j = 0; j < N; j ++) {
      double s = 0.0;
      for (int k = 0; k < N; k ++) {
	s += Q(i, k) * S(k, j);
      }
      
      work(i, j) = s;
    }
  }

  //
  // Compute work Z^T
  //
  for (int i = 0; i < N; i ++) {
    
    for (int j = 0; j < N; j ++) {
      double s = 0.0;
      for (int k = 0; k < N; k ++) {
	s += work(i, k) * Z(j, k);
      }
      
      A(i, j) = s;
    }
  }
}

#endif // specializedeigenproblem_hpp

