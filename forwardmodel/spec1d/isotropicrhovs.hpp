#pragma once
#ifndef isotropicrhovs_hpp
#define isotropicrhovs_hpp

#include "parameterset.hpp"

#include "encodedecode.hpp"

//
// For Love waves, and Rayleigh waves with constant VP/VS ratio
//

template
<
  typename real,
  size_t vpvs_numerator = 17,
  size_t vpvs_denominator = 10
>
class IsotropicRhoVs : public ParameterSet<real, 2> {
public:

  IsotropicRhoVs() :
    ParameterSet<real, 2>()
  {
    prho() = 0.0;
    pVs() = 0.0;
  }

  IsotropicRhoVs(real rho,
		 real Vs) :
    ParameterSet<real, 2>()
  {
    prho() = rho;
    pVs() = Vs;
  }
  
  virtual real rho(real depth) const
  {
    return prho();
  }
  
  virtual real A(real depth) const
  {
    real vpvs = (real)vpvs_numerator/(real)vpvs_denominator;
    
    return prho() * (vpvs * pVs() * vpvs * pVs());
  }
  
  virtual real C(real depth) const
  {
    real vpvs = (real)vpvs_numerator/(real)vpvs_denominator;
    
    return prho() * (vpvs * pVs() * vpvs * pVs());
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
    if (index == 0) {
      return 1.0;
    } else {
      return 0.0;
    }
  }
  
  virtual real dA(size_t index, real depth) const
  {
    real vpvs = (real)vpvs_numerator/(real)vpvs_denominator;
    switch (index) {
    case 0:
      return (vpvs * pVs() * vpvs * pVs());
    case 1:
      return 2.0*vpvs*prho() * pVs();
    default:
      return 0.0;
    }
  }
  
  virtual real dC(size_t index, real depth) const
  {
    real vpvs = (real)vpvs_numerator/(real)vpvs_denominator;
    switch (index) {
    case 0:
      return (vpvs * pVs() * vpvs * pVs());
    case 1:
      return 2.0*vpvs*prho() * pVs();
    default:
      return 0.0;
    }
  }
  
  virtual real dF(size_t index, real depth) const
  {
    real vpvs = (real)vpvs_numerator/(real)vpvs_denominator;
    switch(index) {
    case 0:
      return (vpvs * pVs() * vpvs * pVs()) - 2.0*pVs()*pVs();
    case 1:
      return 2.0*vpvs*prho() * pVs() - 4.0*prho() * pVs();
    default:
      return 0.0;
    }
  }
  
  virtual real dL(size_t index, real depth) const
  {
    switch (index) {
    case 0:
      return pVs()*pVs();
    case 1:
      return 2.0 * prho() * pVs();
    default:
      return 0.0;
    }
  }
  
  virtual real dN(size_t index, real depth) const
  {
    switch (index) {
    case 0:
      return pVs()*pVs();
    case 1:
      return 2.0 * prho() * pVs();
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

  
  real &pVs()
  {
    return this->operator[](1);
  }

  const real &pVs() const
  {
    return this->operator[](1);
  }
};

#endif // isotropicrhovs_hpp
