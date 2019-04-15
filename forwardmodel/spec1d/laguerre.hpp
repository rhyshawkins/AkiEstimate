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
  
