#pragma once
#ifndef symmetricsolve_hpp
#define symmetricsolve_hpp


extern "C" void dsysv_(const char *uplo,
		       int *N,
		       int *NRHS,
		       double *A,
		       int *LDA,
		       int *IPIV,
		       double *B,
		       int *LDB,
		       double *WORK,
		       int *LWORK,
		       int *INFO);

template
<
  typename real
>
bool SymmetricSolve(Spec1DMatrix<real> &A,
		    Spec1DMatrix<real> &B,
		    Spec1DMatrix<real> &work,
		    int *IPIV)
{
  FATAL("Unimplemented");
  return false;
}

template
<>
bool SymmetricSolve<double>(Spec1DMatrix<double> &A,
			    Spec1DMatrix<double> &B,
			    Spec1DMatrix<double> &WORK,
			    int *IPIV)
{
  int N = A.rows();
  int LDA = A.rows();
  int NRHS = B.cols();
  int LDB = B.rows();
  
  int INFO;
  
  int LWORK = -1;
  double tWORK;  
  dsysv_("U", &N, &NRHS, A.data(), &LDA, IPIV, B.data(), &LDB, &tWORK, &LWORK, &INFO);

  if (INFO != 0) {
    FATAL("Failed to get optimal work size");
  }
  
  WORK.resize(tWORK, 1);
  LWORK = WORK.rows();
  
  dsysv_("U", &N, &NRHS, A.data(), &LDA, IPIV, B.data(), &LDB, WORK.data(), &LWORK, &INFO);

  if (INFO != 0) {
    return false;
  }

  return true;
}

#endif // symmetricsolve_hpp
