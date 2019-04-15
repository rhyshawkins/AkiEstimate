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
#ifndef isotropicvsxi_hpp
#define isotropicvsxi_hpp

#include "parameterset.hpp"

template
<
  typename real,
  typename empiricalmodel
>
class IsotropicVsXi : public ParameterSet<real, 2> {
public:

  IsotropicVsXi() :
    ParameterSet<real, 2>()
  {
    pVs() = 0.0;
    pXi() = 0.0;
  }
  
  IsotropicVsXi(real Vs) :
    ParameterSet<real, 2>()
  {
    pVs() = Vs;
    pXi() = 0.0;
  }

  virtual real rho(real depth) const
  {
    real rho, Vp;
    empiricalmodel::compute(pVs(), rho, Vp);
    return rho;
  }
  
  virtual real A(real depth) const
  {
    real rho, Vp;
    empiricalmodel::compute(pVs(), rho, Vp);

    return rho * Vp * Vp;
  }
  
  virtual real C(real depth) const
  {
    real rho, Vp;
    empiricalmodel::compute(pVs(), rho, Vp);

    return rho * Vp * Vp;
  }
  
  virtual real F(real depth) const
  {
    return A(depth) - 2.0 * L(depth);
  }
  
  virtual real L(real depth) const
  {
    real rho, Vp;
    empiricalmodel::compute(pVs(), rho, Vp);

    return rho * pVs() * pVs();
  }
  
  virtual real N(real depth) const
  {
    real rho, Vp;
    empiricalmodel::compute(pVs(), rho, Vp);

    
    return rho * pVs() * pVs() * pXi() * pXi();
  }

  virtual real drho(size_t i, real depth) const
  {
    if (i == 0) {
      real rho, Vp, drho, dvp;
      empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);
      
      return drho;
    } else {
      return 0.0;
    }
  }
  
  virtual real dA(size_t i, real depth) const
  {
    if (i == 0) {
      real rho, Vp, drho, dvp;
      empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);
      
      // A = rho(Vs) * Vp(Vs)^2, want dA/dVs
      return drho * Vp * Vp + rho * 2.0 * Vp * dvp;
    } else {
      return 0.0;
    }
  }
  
  virtual real dC(size_t i, real depth) const
  {
    if (i == 0) {
      real rho, Vp, drho, dvp;
      empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);
      
      // C = rho(Vs) * Vp(Vs)^2, want dA/dVs
      return drho * Vp * Vp + rho * 2.0 * Vp * dvp;
    } else {
      return 0.0;
    }
  }
  
  virtual real dF(size_t i, real depth) const
  {
    if (i == 0) {
      real rho, Vp, drho, dvp;
      empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);
      
      // F = A - 2L
      //   = rho Vp^2 - 2 rho Vs^2
      // 
      // C = rho(Vs) * Vp(Vs)^2, want dA/dVs
      return drho*Vp*Vp + rho*2.0*Vp*dvp - 2.0*(drho*pVs() + rho*2.0)*pVs();
    } else {
      return 0.0;
    }
  }
  
  virtual real dL(size_t i, real depth) const
  {
    if (i == 0) {
      real rho, Vp, drho, dvp;
      empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);
      
      // L = rho(Vs) * Vs^2, want dL/dVs
      return (drho*pVs() + rho*2.0)*pVs();
    } else {
      return 0.0;
    }
  }
  
  virtual real dN(size_t i, real depth) const
  {
    if (i == 0.0) {
      real rho, Vp, drho, dvp;
      empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);
      
      // N = rho(Vs) * (Vs * Xi)^2, want dL/dVs
      return (drho*pVs() + rho*2.0)*pVs()*pXi()*pXi();
    } else {
      real rho, Vp;
      empiricalmodel::compute(pVs(), rho, Vp);

      
      // N = rho(Vs) * (Vs * Xi)^2, want dL/dVs
      return 2.0*rho*pVs()*pVs()*pXi();
    }
  }
  
  real &pVs()
  {
    return this->operator[](0);
  }

  const real &pVs() const
  {
    return this->operator[](0);
  }

  real &pXi()
  {
    return this->operator[](1);
  }

  const real &pXi() const
  {
    return this->operator[](1);
  }
};

#endif // isotropicvsxi_hpp
