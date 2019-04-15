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
#ifndef laguerrequadrature_hpp
#define laguerrequadrature_hpp

#include <array>

#include "laguerre.hpp"
#include "polynomial.hpp"
#include "eigenroots.hpp"

#include "logging.hpp"

//
// This is strictly speaking Laguerre _function_ quadrature, see:
// Jie Shen, "Stable and efficient spectral methods in unbounded domains using Laguerre functions",
// SIAM J. NUMER. ANAL. 38(4), 2000
//
// Nodes are roots of xLn'(x) = x_0 .. x_n
// Weights are 1/((n + 1)(\hat{Ln}(x_i))^2) = w_0 .. w_n
// Cardinal functions are (exp(-x/2)/exp(-x_i/2)) * (-x * L'_{n+1}(x))/((n + 1)L_{n + 1}(x_i)(x - x_i))
// Ln is the nth degree Laguerre polynomial
// \hat{Ln} is the nth degree Laguerre function = Ln exp(-x/2)
//

template
<
  typename real
>
class LaguerreCardinal {
public:

  LaguerreCardinal() :
    scale(1.0),
    P(0, {1.0})
  {
  }
  
  void set(size_t n, size_t i, const real *nodes)
  {
    if (i > n) {
      FATAL("Index out of range");
    }

    if (n == 0) {
      //
      // 0th order
      //
      scale = 1.0;

    } else {

      size_t s = 1;
      for (size_t j = 0; j <= n; j ++) {
	if (j != i) {
	  P *= Polynomial<real>(1, {-nodes[j], 1.0})/(-(real)(s));
	  s ++;
	}
      }

      DynamicLaguerre<real> Ln1(n + 1);

      scale = exp(nodes[i]/2.0) / ((real)(n + 1) * Ln1.value(nodes[i]));
    }

  }

  real value(real x)
  {
    return scale * exp(-x/2.0) * P.value(x);
  }

  real dydx(real x)
  {
    real Pxprime;

    real Px = P.value_and_dydx(x, Pxprime);

    //
    // The difference Pxprime - Px/2 may cause inaccuracies but I think it's ok over
    // the range of the nodes.
    //
    return scale * exp(-x/2.0) * (Pxprime - 0.5*Px);
  }

  real value_and_dydx(real x, real &dydx)
  {
    real Pxprime;

    real Px = P.value_and_dydx(x, Pxprime);

    dydx = scale * exp(-x/2.0) * (Pxprime - 0.5*Px);
    return scale * exp(-x/2.0) * Px;
  }

  void print(FILE *fp)
  {
    fprintf(fp, "%f exp(-x/2) * ", scale);
    P.print(fp);
  }
  

  //
  // Cardinal function is exp(-x/2) * scale * P(x)
  //
  real scale;
  Polynomial<real> P;
};

template
<
  typename real,
  size_t maxorder
>
class LaguerreQuadrature {
public:

  LaguerreQuadrature(size_t order = maxorder) :
    n(order)
  {
    if (n > maxorder) {
      FATAL("Order out of range %d (%d)", (int)n, (int)maxorder);
    }

    switch (n) {
    case 0: // Perhaps shouldn't allow 0th order as the integral is undefined
      nodes[0] = 0.0;
      weights[0] = 1.0;
      cardinal[0].set(0, 0, nodes.data());

      derivative_weights[0][0] = cardinal[0].dydx(nodes[0]);
      break;

    case 1:
      nodes[0] = 0.0;
      nodes[1] = 2.0;

      weights[0] = 1.0/2.0;
      weights[1] = 1.0/(2.0 * exp(-2.0));

      cardinal[0].set(1, 0, nodes.data());
      cardinal[1].set(1, 1, nodes.data());

      derivative_weights[0][0] = cardinal[0].dydx(nodes[0]);
      derivative_weights[0][1] = cardinal[0].dydx(nodes[1]);
      derivative_weights[1][0] = cardinal[1].dydx(nodes[0]);
      derivative_weights[1][1] = cardinal[1].dydx(nodes[1]);
      
      break;

    case 2:

      nodes[0] = 0.0;
      nodes[1] = 3.0 - sqrt(3.0);
      nodes[2] = 3.0 + sqrt(3.0);

      
      weights[0] = 1.0/3.0;
      weights[1] = 1.0/(3.0 * exp(-(3.0 - sqrt(3.0))) * (1.0 - sqrt(3.0))*(1.0 - sqrt(3.0)));
      weights[2] = 1.0/(3.0 * exp(-(3.0 + sqrt(3.0))) * (1.0 + sqrt(3.0))*(1.0 + sqrt(3.0)));

      cardinal[0].set(2, 0, nodes.data());
      cardinal[1].set(2, 1, nodes.data());
      cardinal[2].set(2, 2, nodes.data());

      derivative_weights[0][0] = cardinal[0].dydx(nodes[0]);
      derivative_weights[0][1] = cardinal[0].dydx(nodes[1]);
      derivative_weights[0][2] = cardinal[0].dydx(nodes[2]);
      derivative_weights[1][0] = cardinal[1].dydx(nodes[0]);
      derivative_weights[1][1] = cardinal[1].dydx(nodes[1]);
      derivative_weights[1][2] = cardinal[1].dydx(nodes[2]);
      derivative_weights[2][0] = cardinal[2].dydx(nodes[0]);
      derivative_weights[2][1] = cardinal[2].dydx(nodes[1]);
      derivative_weights[2][2] = cardinal[2].dydx(nodes[2]);
      break;
      
    default:
      {
	//
	// First node always 0
	//
	nodes[0] = 0.0;

	std::vector<real> er;
	int ne;

	if (!eigensolveroots_laguerre(n,
				      er,
				      ne)) {
	  FATAL("Failed to compute eigen values for laguerre nodes");
	}

	if (ne != (int)n) {
	  FATAL("Failed to compute all eigen values for laguerre nodes");
	}
	size_t i = 1;
	for (auto r: er) {
	  nodes[i] = r;
	  i ++;
	}
	
	//
	// Weights
	//
	DynamicLaguerre<real> Ln(n);

	for (size_t i = 0; i <= n; i ++) {
	  real y = Ln.value(nodes[i]);
	  weights[i] = 1.0/((real)(n + 1) * y*y * exp(-nodes[i]));
	}
	
	//
	// Cardinal functions
	//
	for (size_t i = 0; i <= n; i ++) {
	  cardinal[i].set(n, i, nodes.data());
	}
	
	//
	// Derivative weights
	//
	DynamicLaguerre<real> Ln1(n + 1);
	for (size_t i = 0; i <= n; i ++) {
	  for (size_t j = 0; j <= n; j ++) {

	    //
	    // From Shen 3.17
	    //
	    if (j == i) {

	      if (j == 0) {
		derivative_weights[i][j] = -((real)n + 1.0)/2.0;
	      } else {
		derivative_weights[i][j] = 0.0;
	      }

	    } else {
	      derivative_weights[i][j] = (Ln1.value(nodes[j]) * exp(-nodes[j]/2.0))/
		(Ln1.value(nodes[i]) * exp(-nodes[i]/2.0) * (nodes[j] - nodes[i]));
	    }
	  }
	}

      }
      break;
    }
  }

  size_t n;
  std::array<real, maxorder + 1> nodes;
  std::array<real, maxorder + 1> weights;
  std::array<std::array<real, maxorder + 1>, maxorder + 1> derivative_weights;
  std::array<LaguerreCardinal<real>, maxorder + 1> cardinal;
};

#endif // laguerrequadrature_hpp
