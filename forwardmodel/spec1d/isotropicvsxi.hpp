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
