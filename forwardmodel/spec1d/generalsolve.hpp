//
//    Spec1D : A spectral element code for surface wave dispersion of Love
//    and Rayleigh waves. See
//
//      R Hawkins, "A spectral element method for surface wave dispersion and adjoints",
//      Geophysical Journal International, 2018, 215:1, 267 - 302
//      https://doi.org/10.1093/gji/ggy277
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//


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
