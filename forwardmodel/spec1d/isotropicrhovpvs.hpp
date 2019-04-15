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
#ifndef isotropicrhovpvs_hpp
#define isotropicrhovpvs_hpp

#include "parameterset.hpp"

template
<
  typename real
>
class IsotropicRhoVpVs : public ParameterSet<real, 3> {
public:

  IsotropicRhoVpVs() :
    ParameterSet<real, 3>()
  {
    prho() = 0.0;
    pVp() = 0.0;
    pVs() = 0.0;
  }

  IsotropicRhoVpVs(real rho,
		   real Vp,
		   real Vs) :
    ParameterSet<real, 3>()
  {
    prho() = rho;
    pVp() = Vp;
    pVs() = Vs;
  }
  
  virtual real rho(real depth) const
  {
    return prho();
  }
  
  virtual real A(real depth) const
  {
    return prho() * pVp() * pVp();
  }
  
  virtual real C(real depth) const
  {
    return prho() * pVp() * pVp();
  }
  
  virtual real F(real depth) const
  {
    return A(depth) - 2.0 * L(depth);
  }
  
  virtual real L(real depth) const
  {
    return prho() * pVs() * pVs();
  }
  
  virtual real N(real depth) const
  {
    return prho() * pVs() * pVs();
  }

  virtual real drho(size_t index, real depth) const
  {
    switch (index) {
    case 0:
      return 1.0;
    default:
      return 0.0;
    }
  }
  
  virtual real dA(size_t index, real depth) const
  {
    // A = rho Vp^2
    switch (index) {
    case 0:
      return pVp()*pVp();
    case 1:
      return 2.0*prho()*pVp();
    case 2:
    default:
      return 0.0;
    }
  }
  
  virtual real dC(size_t index, real depth) const
  {
    // C = rho Vp^2
    switch (index) {
    case 0:
      return pVp()*pVp();
    case 1:
      return 2.0*prho()*pVp();
    case 2:
    default:
      return 0.0;
    }
  }
  
  virtual real dF(size_t index, real depth) const
  {
    // F = A - 2L
    //   = rhoVp^2 - 2rhoVs^2
    switch (index) {
    case 0:
      return pVp()*pVp() - 2.0*pVs()*pVs();
    case 1:
      return 2.0*prho()*pVp();
    case 2:
      return -4.0*prho()*pVs();
    default:
      return 0.0;
    }
  }
  
  virtual real dL(size_t index, real depth) const
  {
    // L = rho Vs^2
    switch (index) {
    case 0:
      return pVs()*pVs();
    case 2:
      return 2.0*prho()*pVs();
    case 1:
    default:
      return 0.0;
    }
  }
  
  virtual real dN(size_t index, real depth) const
  {
    // N = rho Vs^2
    switch (index) {
    case 0:
      return pVs()*pVs();
    case 2:
      return 2.0*prho()*pVs();
    case 1:
    default:
      return 0.0;
    }
  }

private:

  real &prho()
  {
    return this->operator[](0);
  }
  
  const real &prho() const
  {
    return this->operator[](0);
  }

  real &pVp()
  {
    return this->operator[](1);
  }
  
  const real &pVp() const
  {
    return this->operator[](1);
  }

  real &pVs()
  {
    return this->operator[](2);
  }
  
  const real &pVs() const
  {
    return this->operator[](2);
  }
};

#endif // isotropicrhovpvs_hpp
