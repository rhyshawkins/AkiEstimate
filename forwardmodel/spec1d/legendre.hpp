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

    
