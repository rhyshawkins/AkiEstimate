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
