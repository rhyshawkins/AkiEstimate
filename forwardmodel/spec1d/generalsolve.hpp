#pragma once
#ifndef generalsolve_hpp
#define generalsolve_hpp

extern "C" void dgesv_(int *N,
		       int *NRHS,
		       double *A,
		       int *LDA,
		       int *IPIV,
		       double *B,
		       int *LDB,
		       int *INFO);

template
<
  typename real
>
bool GeneralSolve(Spec1DMatrix<real> &A,
		  Spec1DMatrix<real> &B,
		  int *IPIV)
{
  FATAL("Unimplemented");
  return false;
}

template
<>
bool GeneralSolve<double>(Spec1DMatrix<double> &A,
			  Spec1DMatrix<double> &B,
			  int *IPIV)
{
  int N = A.rows();
  int LDA = A.rows();
  int NRHS = B.cols();
  int LDB = B.rows();
  
  int INFO;
  
  
  dgesv_(&N, &NRHS, A.data(), &LDA, IPIV, B.data(), &LDB, &INFO);

  if (INFO != 0) {
    return false;
  }

  return true;
}

#endif // generalsolve_hpp
