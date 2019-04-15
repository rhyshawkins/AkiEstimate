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
#ifndef regression_hpp
#define regression_hpp

#include "parameterset.hpp"

#include "encodedecode.hpp"

template
<
  typename real
>
class Regression : public ParameterSet<real, 1> {
public:

  Regression() :
    ParameterSet<real, 1>()
  {
    y() = 0.0;
  }

  Regression(real _y) :
    ParameterSet<real, 1>()
  {
    y() = _y;
  }
  
  virtual real rho(real depth) const
  {
    return y();
  }
  
  virtual real A(real depth) const
  {
    return y();
  }
  
  virtual real C(real depth) const
  {
    return y();
  }
  
  virtual real F(real depth) const
  {
    return y();
  }
  
  virtual real L(real depth) const
  {
    return y();
  }
  
  virtual real N(real depth) const
  {
    return y();
  }

  virtual real drho(size_t index, real depth) const
  {
    return 1.0;
  }
  
  virtual real dA(size_t index, real depth) const
  {
    return 1.0;
  }
  
  virtual real dC(size_t index, real depth) const
  {
    return 1.0;
  }
  
  virtual real dF(size_t index, real depth) const
  {
    return 1.0;
  }
  
  virtual real dL(size_t index, real depth) const
  {
    return 1.0;
  }
  
  virtual real dN(size_t index, real depth) const
  {
    return 1.0;
  }

  static const char *NAME()
  {
    return "Regression";
  }
  
private:

  real &y()
  {
    return this->operator[](0);
  }

  const real &y() const
  {
    return this->operator[](0);
  }
  
};

#endif // regression_hpp
