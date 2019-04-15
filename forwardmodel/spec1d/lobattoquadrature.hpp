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
#ifndef lobattoquadrature_hpp
#define lobattoquadrature_hpp

#include <array>

#include "legendre.hpp"
#include "polynomial.hpp"
#include "eigenroots.hpp"

#include "logging.hpp"

template
<
  typename real,
  size_t maxorder
>
class LobattoQuadrature {
public:

  LobattoQuadrature(size_t order = maxorder) :
    n(order)
  {
    if (n > maxorder) {
      FATAL("Order out of range %d (%d)", (int)n, (int)maxorder);
    }


    switch (n) {
    case 0:
      nodes[0] = 0.0;
      weights[0] = 2.0;
      cardinal[0] = Polynomial<real>(0, {1.0});

      derivative_weights[0][0] = 0.0;
      break;

    case 1:
      nodes[0] = -1.0;
      nodes[1] = 1.0;
      
      weights[0] = 1.0;
      weights[1] = 1.0;
      
      cardinal[0] = Polynomial<real>(1, {0.5, -0.5});
      cardinal[1] = Polynomial<real>(1, {0.5, 0.5});

      derivative_weights[0][0] = -0.5;
      derivative_weights[0][1] = -0.5;
      derivative_weights[1][0] = 0.5;
      derivative_weights[1][1] = 0.5;
      
      break;
      
    case 2:
      nodes[0] = -1.0;
      nodes[1] = 0.0;
      nodes[2] = 1.0;
      
      weights[0] = 1.0/3.0;
      weights[1] = 4.0/3.0;
      weights[2] = 1.0/3.0;
      
      cardinal[0] = Polynomial<real>(2, {0.0, -0.5, 0.5});
      cardinal[1] = Polynomial<real>(2, {1.0, 0.0, -1.0});
      cardinal[2] = Polynomial<real>(2, {0.0, 0.5, 0.5});

      derivative_weights[0][0] = -1.5;
      derivative_weights[0][1] = -0.5;
      derivative_weights[0][2] = 0.5;

      derivative_weights[1][0] = 2.0;
      derivative_weights[1][1] = 0.0;
      derivative_weights[1][2] = -2.0;

      derivative_weights[2][0] = -0.5;
      derivative_weights[2][1] = 0.5;
      derivative_weights[2][2] = 1.5;
      break;

    default:
      {
	//
	// End nodes
	//
	nodes[0] = -1.0;
	nodes[n] = 1.0;

	std::vector<real> er;
	int ne;

	if (!eigensolveroots_lobatto(n,
				     er,
				     ne)) {
	  FATAL("Failed to compute eigen values for lobatto nodes");
	}

	if (ne != (int)(n - 1)) {
	  FATAL("Failed to compute all eigen values for lobatto nodes");
	}
	size_t i = 1;
	for (auto r: er) {
	  nodes[i] = r;
	  i ++;
	}
	
	//
	// End weights
	//
	weights[0] = 2.0/((real)n * (real)(n + 1));
	weights[n] = weights[0];

	
	DynamicLegendre<real> Pn(n);
	
	for (size_t i = 1; i < n; i ++) {
	  real y = Pn.value(nodes[i]);
	  weights[i] = 2.0/(y*y * (real)n * (real)(n + 1));
	}
	
	//
	// Cardinal functions
	//
	for (size_t i = 0; i <= n; i ++) {
	  cardinal[i] = Polynomial<real>(n, i, nodes.data());
	}

	//
	// Derivative weights
	//
	for (size_t i = 0; i <= n; i ++) {
	  Polynomial<real> lprime = cardinal[i].derivative();

	  for (size_t j = 0; j <= n; j ++) {

	    derivative_weights[i][j] = lprime.value(nodes[j]);

	  }
	}
      }
      break;
    }
  }

  size_t n;
  std::array<real, maxorder + 1> nodes;
  std::array<real, maxorder + 1> weights;

  // Derivative weights[i][j] are cardinal_i'(nodes[j])
  std::array<std::array<real, maxorder + 1>, maxorder + 1> derivative_weights;
  
  std::array<Polynomial<real>, maxorder + 1> cardinal;

};

#endif // lobattoquadrature_hpp
