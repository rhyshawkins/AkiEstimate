#pragma once
#ifndef lobatto_hpp
#define lobatto_hpp

#include "polynomial.hpp"
#include "legendre.hpp"
template
<
  typename real,
  size_t order
>
class Lobatto : public Polynomial {
public:

  Lobatto() :
    Polynomial(Legendre<real, order + 1>.derivative())
  {
  }

};

#endif // lobatto_hpp
    
