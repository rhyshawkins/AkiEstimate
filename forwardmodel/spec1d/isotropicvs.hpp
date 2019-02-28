#pragma once
#ifndef isotropicvs_hpp
#define isotropicvs_hpp

#include "parameterset.hpp"

template
<
  typename real,
  typename empiricalmodel
>
class IsotropicVs : public ParameterSet<real, 1> {
public:

  IsotropicVs() :
    ParameterSet<real, 1>()
  {
    pVs() = 0.0;
  }
  
  IsotropicVs(real Vs) :
    ParameterSet<real, 1>()
  {
    pVs() = Vs;
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

    return rho * pVs() * pVs();
  }

  virtual real drho(size_t i, real depth) const
  {
    real rho, Vp, drho, dvp;
    empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);

    return drho;
  }
  
  virtual real dA(size_t i, real depth) const
  {
    real rho, Vp, drho, dvp;
    empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);

    // A = rho(Vs) * Vp(Vs)^2, want dA/dVs
    return drho * Vp * Vp + rho * 2.0 * Vp * dvp;
  }
  
  virtual real dC(size_t i, real depth) const
  {
    real rho, Vp, drho, dvp;
    empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);

    // C = rho(Vs) * Vp(Vs)^2, want dA/dVs
    return drho * Vp * Vp + rho * 2.0 * Vp * dvp;
  }
  
  virtual real dF(size_t i, real depth) const
  {
    real rho, Vp, drho, dvp;
    empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);

    // F = A - 2L
    //   = rho Vp^2 - 2 rho Vs^2
    // 
    // C = rho(Vs) * Vp(Vs)^2, want dA/dVs
    return drho*Vp*Vp + rho*2.0*Vp*dvp - 2.0*(drho*pVs() + rho*2.0)*pVs();
  }
  
  virtual real dL(size_t i, real depth) const
  {
    real rho, Vp, drho, dvp;
    empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);

    // L = rho(Vs) * Vs^2, want dL/dVs
    return (drho*pVs() + rho*2.0)*pVs();
  }
  
  virtual real dN(size_t i, real depth) const
  {
    real rho, Vp, drho, dvp;
    empiricalmodel::compute_gradient(pVs(), rho, Vp, drho, dvp);

    // N = rho(Vs) * Vs^2, want dL/dVs
    return (drho*pVs() + rho*2.0)*pVs();
  }
  
  real &pVs()
  {
    return this->operator[](0);
  }

  const real &pVs() const
  {
    return this->operator[](0);
  }

  static const char *NAME()
  {
    return "IsotropicVs";
  }
};

#endif // isotropicvs_hpp
