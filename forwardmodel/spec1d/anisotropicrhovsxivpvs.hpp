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
#ifndef anisotropicrhovsxivpvs_hpp
#define anisotropicrhovsxivpvs_hpp

//
// Parameterisation of
//  rho
//  vs = sqrt((2L + N)/(3*rho))
//  xi = N/L
//  vpvs
//

//
// Mappings
//
// rho -> rho
// A   -> rho * (vs * vpvs)^2
// C   -> same
// F   -> A - 2L
// L   -> 3 rho vs^2 / (2 + xi)
// N   -> 3 rho vs^2 / (2/xi + 1)
//
#include "parameterset.hpp"

template
<
  typename real
>
class AnisotropicRhoVsXiVpVs : public ParameterSet<real, 4> {
public:

  AnisotropicRhoVsXiVpVs() :
    ParameterSet<real, 4>()
  {
    prho() = 0.0;
    pVs() = 0.0;
    pXi() = 0.0;
    pVpVs() = 0.0;
    
  }

  AnisotropicRhoVsXiVpVs(real rho,
			 real Vs,
			 real Xi,
			 real VpVs) :
    ParameterSet<real, 4>()
  {
    prho() = rho;
    pVs() = Vs;
    pXi() = Xi;
    pVpVs() = VpVs;
  }
  
  virtual real rho(real depth) const
  {
    return prho();
  }
  
  virtual real A(real depth) const
  {
    real vp = pVs() * pVpVs();
    return prho() * vp * vp;
  }
  
  virtual real C(real depth) const
  {
    real vp = pVs() * pVpVs();
    return prho() * vp * vp;
  }
  
  virtual real F(real depth) const
  {
    return A(depth) - 2.0 * L(depth);
  }
  
  virtual real L(real depth) const
  {
    return 3.0 * prho() * pVs() * pVs()/(2.0 + pXi());
  }
  
  virtual real N(real depth) const
  {
    return 3.0 * prho() * pVs() * pVs()/(2.0/pXi() + 1.0);
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
    // A = vs^2 vpvs^2 rho
    switch (index) {
    case 0: // rho
      return pVs()*pVs()*pVpVs()*pVpVs();
    case 1: // vs
      return 2.0*pVs()*prho()*pVpVs()*pVpVs();
    case 2: // xi
      return 0.0;
    case 3: // vpvs
      return 2.0*pVs()*pVs()*prho()*pVpVs();
    default:
      return 0.0;
    }
  }
  
  virtual real dC(size_t index, real depth) const
  {
    // C = vs^2 vpvs^2 rho
    switch (index) {
    case 0: // rho
      return pVs()*pVs()*pVpVs()*pVpVs();
    case 1: // vs
      return 2.0*pVs()*prho()*pVpVs()*pVpVs();
    case 2: // xi
      return 0.0;
    case 3: // vpvs
      return 2.0*pVs()*pVs()*prho()*pVpVs();
    default:
      return 0.0;
    }
  }
  
  virtual real dF(size_t index, real depth) const
  {
    // F = A - 2L
    return dA(index, depth) - 2.0*dL(index, depth);
  }
  
  virtual real dL(size_t index, real depth) const
  {
    // L = 3 rho Vs^2/(2 + xi)
    switch (index) {
    case 0: // rho
      return 3.0 * pVs()*pVs()/(2.0 + pXi());
    case 1: // vs
      return 6.0 * prho()* pVs()/(2.0 + pXi());
    case 2: // xi
      return -3.0 * prho() * pVs() * pVs() / ((2.0 + pXi()) * (2.0 + pXi()));
    case 3: // vpvs
      return 0.0;
    default:
      return 0.0;
    }
  }
  
  virtual real dN(size_t index, real depth) const
  {
    // N = 3 rho vs^2/(2/xi + 1)
    switch (index) {
    case 0: // rho
      return 3.0 * pVs()*pVs()/(2.0/pXi() + 1);
    case 1: // vs
      return 6.0 * prho()* pVs()/(2.0/pXi() + 1);
    case 2: // xi
      return 6.0 * prho() * pVs() * pVs() / ((2.0 + pXi()) * (2.0 + pXi()));
    case 3: // vpvs
      return 0.0;
    default:
      return 0.0;
    }
  }

  static const char *NAME()
  {
    return "AnisotropicRhoVsXiVpVs";
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

  real &pVs()
  {
    return this->operator[](1);
  }
  
  const real &pVs() const
  {
    return this->operator[](1);
  }

  real &pXi()
  {
    return this->operator[](2);
  }
  
  const real &pXi() const
  {
    return this->operator[](2);
  }

  real &pVpVs()
  {
    return this->operator[](3);
  }
  
  const real &pVpVs() const
  {
    return this->operator[](3);
  }
  
};


#endif // anisotropicrhovsxivpvs_hpp
