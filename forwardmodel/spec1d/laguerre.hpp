#pragma once
#ifndef laguerre_hpp
#define laguerre_hpp

#include "polynomial.hpp"

template
<
  typename real,
  size_t n
>
class Laguerre : public Polynomial<real> {
public:

  Laguerre() :
    Polynomial<real>((Polynomial<real>(1, {(real)n*2.0 - 1.0, -1.0}) * Laguerre<real, n - 1>() -
		      (real)(n - 1)*Laguerre<real, n - 2>())/(real)(n))
  {
  }
};

template
<
  typename real
>
class Laguerre<real, 0> : public Polynomial<real> {
public:

  Laguerre() :
    Polynomial<real>(0, {1.0})
  {
  }
  
};

template
<
  typename real
>
class Laguerre<real, 1> : public Polynomial<real> {
public:

  Laguerre() :
    Polynomial<real>(1, {1.0, -1.0})
  {
  }
  
};

template
<
  typename real
>
class DynamicLaguerre : public Polynomial<real> {
public:

  static constexpr size_t CACHE = 100;
  static std::array<Polynomial<real>, CACHE> cache;
  
  DynamicLaguerre(size_t n)
  {
    switch (n) {
    case 0:
      *this = Polynomial<real>(0, {1.0});
      break;

    case 1:
      *this = Polynomial<real>(1, {1.0, -1.0});
      break;

    default:
      if (n < CACHE) {
	//
	// Memoization of polynomials to speed up computation for large n
	//
	if (cache[n].order() == 0) {
	  cache[n] = (Polynomial<real>(1, {(real)n*2.0 - 1.0, -1.0}) * DynamicLaguerre(n - 1) -
	       (real)(n - 1)*DynamicLaguerre(n - 2))/(real)(n);
	}
	*this = cache[n];
      } else {
	*this = (Polynomial<real>(1, {(real)n*2.0 - 1.0, -1.0}) * DynamicLaguerre(n - 1) -
		 (real)(n - 1)*DynamicLaguerre(n - 2))/(real)(n);
      }

    }
  }

  DynamicLaguerre& operator=(const Polynomial<real> &rhs)
  {
    Polynomial<real>::operator=(rhs);
    return *this;
  }
    
  DynamicLaguerre& operator=(const DynamicLaguerre &rhs)
  {
    Polynomial<real>::operator=((const Polynomial<real> &)rhs);
    return *this;
  }
};

template <typename real> std::array<Polynomial<real>, DynamicLaguerre<real>::CACHE> DynamicLaguerre<real>::cache;

#endif // laguerre_hpp
  
