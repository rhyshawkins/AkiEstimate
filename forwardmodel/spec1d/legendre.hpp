#pragma once
#ifndef legendre_hpp
#define legendre_hpp

#include "polynomial.hpp"

template
<
  typename real,
  size_t n
>
class Legendre : public Polynomial<real> {
public:

  Legendre() :
    Polynomial<real>((2.0*n - 1.0)/n * Legendre<real, 1>() * Legendre<real, n - 1>() -
		     (n - 1.0)/n * Legendre<real, n - 2>())
  {
  }

};

template
<
  typename real
>
class Legendre<real, 0> : public Polynomial<real> {
public:

  Legendre() :
    Polynomial<real>(0, {1.0})
  {
  }
  
};
  
template
<
  typename real
>
class Legendre<real, 1> : public Polynomial<real> {
public:

  Legendre() :
    Polynomial<real>(1, {0.0, 1.0})
  {
  }
  
};

template
<
  typename real
>
class DynamicLegendre : public Polynomial<real> {
public:

  DynamicLegendre(size_t n)
  {
    switch(n) {
    case 0:
      *this = Polynomial<real>(0, {1.0});
      break;

    case 1:
      *this = Polynomial<real>(1, {0.0, 1.0});
      break;

    default:
      *this = (2.0*n - 1.0)/n * DynamicLegendre<real>(1) * DynamicLegendre<real>(n - 1) -
	(n - 1.0)/n * DynamicLegendre<real>(n - 2);
    }
  }

  DynamicLegendre& operator=(const Polynomial<real> &rhs)
  {
    Polynomial<real>::operator=(rhs);
    return *this;
  }
    
  DynamicLegendre& operator=(const DynamicLegendre &rhs)
  {
    Polynomial<real>::operator=((const Polynomial<real> &)rhs);
    return *this;
  }
};

#endif // legendre_hpp

    
