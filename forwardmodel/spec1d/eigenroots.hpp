#pragma once
#ifndef eigenroots_hpp
#define eigenroots_hpp

#include <vector>
#include <algorithm>

#include "polynomial.hpp"
#include "spec1dmatrix.hpp"
#include "generalisedeigenproblem.hpp"

//
// Solve roots of generic polynomial using eigen values of companion form matrix
//
template
<
  typename real
>
bool
eigensolveroots(const Polynomial<real> &poly,
		std::vector<real> &eig_real,
		std::vector<real> &eig_imag,
		int &nreal)
{
  Spec1DMatrix<real> A;
  Spec1DMatrix<real> I;
  Spec1DMatrix<real> work, eu, ev, lambda;

  int N = poly.order();
  long double denom = poly[N];

  A.resize(N, N);
  A.setZero();

  I.setIdentity(N);

  for (int i = 0; i < N; i ++) {
    A(0, i) = -poly[N - i - 1]/denom;

    if (i < (N - 1)) {
      A(i + 1, i) = 1.0;
    }
  }

  for (int j = 0; j < N; j ++) {
    for (int i = 0; i < N; i ++) {

      printf("%.20lg ", (double)A(j, i));

    }
    printf("\n");
  }

  if (!GEP(A, I, work, eu, ev, lambda)) {
    ERROR("Failed to compute eigen values");
    return false;
  }

  nreal = 0;
  int nimag = 0;

  eig_real.resize(N);
  eig_imag.resize(N);
  
  for (int i = 0; i < N; i ++) {
    if (lambda(i, 1) == 0.0) {
      eig_imag[nreal] = 0.0;
      eig_real[nreal] = lambda(i, 0)/lambda(i, 2);

      nreal ++;

    } else {
      eig_imag[N - 1 - nimag] = lambda(i, 1)/lambda(i, 2);
      eig_real[N - 1 - nimag] = lambda(i, 0)/lambda(i, 2);

      nimag ++;
    }      
  }

  std::sort(eig_real.begin(), eig_real.begin() + nreal);
  
  return true;
}

//
// Solve for internal roots deriviative of Legendre function for a given order for Gauss-Legendre-Lobatto
// nodes
//
template
<
  typename real
>
bool eigensolveroots_lobatto(size_t order,
			     std::vector<real> &eig_real,
			     int &nreal)
{
  switch (order) {
  case 0:
    eig_real.push_back(0.0);
    nreal = 1;
    break;

  case 1:
    eig_real.push_back(-1.0);
    eig_real.push_back(1.0);
    nreal = 2;
    break;

  case 2:
    eig_real.push_back(-1.0);
    eig_real.push_back(0.0);
    eig_real.push_back(1.0);
    nreal = 3;
    break;

  default:
    {
      Spec1DMatrix<double> A;
      Spec1DMatrix<double> I;
      Spec1DMatrix<double> work, eu, ev, lambda;
      
      size_t N = order - 1;
      A.resize(N, N);
      A.setZero();

      I.setIdentity(N);
      
      for (size_t i = 0; i < (N - 1); i ++) {

	if (i == 0) {
	  A(i, i + 1) = 1.0/sqrt(5.0);
	  A(i + 1, i) = 1.0/sqrt(5.0);
	} else {
	  long double ii = 1.0 + i;
	  A(i, i + 1) = 2.0 * sqrt(ii * (ii + 1.0) * (ii + 1.0) * (ii + 2.0)/((2.0*ii + 2.0) * (2.0*ii + 2.0) - 1.0))/(2.0*ii + 2.0);
	  A(i + 1, i) = A(i, i + 1);
	}

      }

      // for (size_t j = 0; j < N; j ++) {
      // 	for (size_t i = 0; i < N; i ++) {

      // 	  printf("%7.3f ", (double)A(j, i));
	  
      // 	}
      // 	printf("\n");
      // }
      if (!GEP(A, I, work, eu, ev, lambda)) {
	ERROR("Failed to compute eigen values");
	return false;
      }
      
      nreal = 0;

      eig_real.resize(N);
  
      for (size_t i = 0; i < N; i ++) {

	if (lambda(i, 1) == 0.0) {
	  eig_real[nreal] = lambda(i, 0)/lambda(i, 2);
	  
	  nreal ++;
	}
      }
      
      std::sort(eig_real.begin(), eig_real.begin() + nreal);

    }
  }

  return true;
}

//
// Solve for roots of the derivative of Laguerre function for a given order for Gauss-Legendre-Laguerre
// nodes
//
template
<
  typename real
>
bool eigensolveroots_laguerre(size_t order,
			      std::vector<real> &eig_real,
			      int &nreal)
{
  switch (order) {
  case 0:
    nreal = 0;
    break;

  case 1:
    eig_real.push_back(2.0);
    nreal = 1;
    break;

  default:
    {
      Spec1DMatrix<real> A;
      Spec1DMatrix<real> I;
      Spec1DMatrix<real> work, eu, ev, lambda;
      
      size_t N = order;
      A.resize(N, N);
      A.setZero();
      I.setIdentity(N);
      
      for (size_t i = 0; i < N; i ++) {

	A(i, i) = 2.0 * (i + 1);

	if (i < (N - 1)) {

	  A(i, i + 1) = sqrt((1.0 + i) * (2.0 + i));
	  A(i + 1, i) = A(i, i + 1);
	  
	}

      }

      if (!GEP(A, I, work, eu, ev, lambda)) {
	ERROR("Failed to compute eigen values");
	return false;
      }

      nreal = 0;

      eig_real.resize(N);

      for (size_t i = 0; i < N; i ++) {
	if (lambda(i, 1) == 0.0) {
	  
	  eig_real[nreal] = lambda(i, 0)/lambda(i, 2);
	
	  nreal ++;

	}
	  
      }
      
      std::sort(eig_real.begin(), eig_real.begin() + nreal);

    }
  }

  return true;
}


#endif // eigenroots_hpp
