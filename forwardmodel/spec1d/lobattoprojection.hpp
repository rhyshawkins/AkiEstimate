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
#ifndef lobattoprojection_hpp
#define lobattoprojection_hpp

#include <array>

#include "lobattoquadrature.hpp"

template
<
  typename real,
  size_t maxorder
>
class LobattoProjection {
public:

  LobattoProjection(size_t order = maxorder)
  {
    //
    // 0th order is straight forward
    //
    for (size_t j = 0; j <= order; j ++) {
      for (size_t k = 0; k <= order; k ++) {
	if (j == 0) {
	  weight[0][j][k] = 1.0;
	} else {
	  weight[0][j][k] = 0.0;
	}
      }
    }

    LobattoQuadrature<real, maxorder> LqHigh(order);
    
    for (size_t i = 1; i < order; i ++) {

      //
      // The weight is the value of the cardinal function at the given higher order node
      // position
      //
      LobattoQuadrature<real, maxorder> LqLow(i);

      for (size_t j = 0; j <= i; j ++) {

	for (size_t k = 0; k <= order; k ++) {

	  weight[i][j][k] = LqLow.cardinal[j].value(LqHigh.nodes[k]);
	  // printf("%d %d %d: %f %f: ", (int)i, (int)j, (int)k, weight[i][j][k], LqHigh.nodes[k]);
	  // LqLow.cardinal[j].print(stdout);
	}
      }

      
    }

    //
    // order -> order is also a direct translation
    //
    for (size_t j = 0; j <= order; j ++) {
      for (size_t k = 0; k <= order; k ++) {
	if (j == k) {
	  weight[order][j][k] = 1.0;
	} else {
	  weight[order][j][k] = 0.0;
	}
      }
    }
  }

  
  //
  // This 3D array is the set of weights to project lower orders to the given
  // template parameter order. weight[i][j][k] is the weight for order
  // i of the jth node of order i to the kth node of the order.
  //
  // With the order template parameter set to 2, ie quadratic, the array
  // will be
  //
  // {
  //   {
  //     { 1, 0, 0 }, { 1, 0, 0 }, { 1, 0, 0 }     0th order to 2nd
  //   }, {
  //     { 1, 0, 0 }, { 0.5, 0.5, 0 }, { 0, 1, 0 } 1st order to 2nd
  //   }, {
  //     { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }     2nd order to 2nd
  //   }
  // }
  //

  
      
  std::array<std::array<std::array<real, maxorder + 1>, maxorder + 1>, maxorder + 1> weight; 
};

template
<
  typename real,
  size_t maxorder
>
class LobattoUpProjection {
public:

  LobattoUpProjection()
  {
  }
  
  LobattoUpProjection(size_t order)
  {
    init(order);
  }

  void init(size_t order)
  {
    if ((order + 1) > maxorder) {
      FATAL("Order out of range");
    }

    LobattoQuadrature<real, maxorder> LqLow(order);
    LobattoQuadrature<real, maxorder> LqHigh(order + 1);
    
    //
    // The weight is the value of the cardinal function at the given higher order node
    // position
    //
    
    for (size_t j = 0; j <= order; j ++) {
      for (size_t k = 0; k <= (order + 1); k ++) {
	weight[j][k] = LqLow.cardinal[j].value(LqHigh.nodes[k]);
      }
    }
  }

  std::array<std::array<real, maxorder + 1>, maxorder + 1> weight;
};

template
<
  typename real,
  size_t maxorder
>
class LobattoDownProjection {
public:

  LobattoDownProjection()
  {
  }
  
  LobattoDownProjection(size_t order)
  {
    init(order);
  }

  void init(size_t order)
  {
    if (order < 1 || order > maxorder) {
      FATAL("Order out of range");
    }

    LobattoQuadrature<real, maxorder> LqLow(order - 1);
    LobattoQuadrature<real, maxorder> LqHigh(order);
    
    //
    // The weight is the value of the cardinal function at the given higher order node
    // position
    //
    
    for (size_t j = 0; j <= order; j ++) {
      for (size_t k = 0; k < order; k ++) {
	weight[j][k] = LqHigh.cardinal[j].value(LqLow.nodes[k]);
      }
    }
  }

  std::array<std::array<real, maxorder + 1>, maxorder + 1> weight;
};

template
<
  typename real,
  size_t maxorder
>
class LobattoUpProjector {
public:

  LobattoUpProjector()
  {
    for (size_t i = 0; i < maxorder; i ++) {
      P[i].init(i);
    }
  }

  void upproject(size_t order, const real *source, real *dest)
  {
    if (order >= maxorder) {
      FATAL("Order out of range");
    }

    for (size_t i = 0; i <= (order + 1); i ++) {

      dest[i] = 0.0;

      for (size_t j = 0; j <= order; j ++) {

	dest[i] += source[j] * P[order].weight[j][i];

      }
    }
  }

  std::array<LobattoUpProjection<real, maxorder>, maxorder> P;
};

template
<
  typename real,
  size_t maxorder
>
class LobattoDownProjector {
public:

  LobattoDownProjector()
  {
    for (size_t i = 0; i < maxorder; i ++) {
      P[i].init(i + 1);
    }
  }

  void downproject(size_t order, const real *source, real *dest)
  {
    if (order == 0 || order > maxorder) {
      FATAL("Order out of range");
    }

    for (size_t i = 0; i < order; i ++) {
      dest[i] = 0.0;

      for (size_t j = 0; j <= order; j ++) {

	dest[i] += source[j] * P[order - 1].weight[j][i];

      }
    }
  }

  
  std::array<LobattoDownProjection<real, maxorder>, maxorder> P;
};

  

#endif // lobattoprojection_hpp
