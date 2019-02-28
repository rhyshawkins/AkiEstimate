#pragma once
#ifndef polynomial_hpp
#define polynomial_hpp

#include <memory>
#include <cmath>

#include "logging.hpp"

template
<
  typename real,
  size_t DEFAULT_SIZE = 16
>
class Polynomial {
public:

  Polynomial() :
    size(DEFAULT_SIZE),
    N(0),
    coeff(new real[DEFAULT_SIZE])
  {
    coeff[0] = 0.0;
  }

  Polynomial(size_t order, const real *coefficients) :
    size(std::max<size_t>(order + 1, DEFAULT_SIZE)),
    N(order),
    coeff(new real[std::max<size_t>(order + 1, DEFAULT_SIZE)])
  {
    for (size_t i = 0; i <= order + 1; i ++) {
      coeff[i] = coefficients[i];
    }
  }

  Polynomial(size_t order, const std::initializer_list<real> &coefficients) :
    size(std::max<size_t>(order + 1, DEFAULT_SIZE)),
    N(order),
    coeff(new real[std::max<size_t>(order + 1, DEFAULT_SIZE)])
  {
    if ((order + 1) != coefficients.size()) {
      FATAL("Initializer size mismatch");
    } else {
      size_t i = 0;
      for (auto c: coefficients) {
	coeff[i] = c;
	i ++;
      }
    }
  }

  //
  // Construct polynomial as a Lagrange cardinal polynomial which is a polynomial
  // that is 1 at x[cardinal] and 0 at other nodal points. The x points are
  // generally Legendre, Lobatto, or Laguerre nodes.
  //
  Polynomial(size_t order,
	     size_t cardinal,
	     const real *x) :
    size(std::max<size_t>(order + 1, DEFAULT_SIZE)),
    N(0),
    coeff(new real[std::max<size_t>(order + 1, DEFAULT_SIZE)])
  {
    if (cardinal > order) {
      FATAL("Cardinal index out of range");
    } else {

      coeff[0] = 1.0;
      for (size_t i = 1; i <= order; i ++) {
	coeff[i] = 0.0;
      }

      for (size_t k = 0; k <= order; k ++) {

	if (k != cardinal) {

	  real denom = x[cardinal] - x[k];
	  (*this) *= Polynomial(1, {-x[k]/denom, 1.0/denom});
	}
      }
    }
  }
  
  Polynomial(const Polynomial &rhs) :
    size(std::max<size_t>(rhs.N + 1, DEFAULT_SIZE)),
    N(rhs.N),
    coeff(new real[std::max<size_t>(rhs.N + 1, DEFAULT_SIZE)])
  {
    for (size_t i = 0; i <= N; i ++) {
      coeff[i] = rhs.coeff[i];
    }
  }

  Polynomial derivative()
  {
    Polynomial p(*this);
    
    for (size_t i = 0; i < p.N; i ++) {
      p.coeff[i] = p.coeff[i + 1] * (real)(i + 1);
    }
    p.N --;

    return p;
  }

  real value(real x)
  {
    real accum = coeff[N];
    for (int i = (int)(N - 1); i >= 0; i --) {
      accum = accum * x + coeff[i];
    }

    return accum;
  }

  real dydx(real x)
  {
    real accum = coeff[N] * (real)N;
    for (int i = (int)(N - 1); i > 0; i --) {
      accum = accum * x + coeff[i] * (real)i;
    }

    return accum;
  }

  real value_and_dydx(real x, real &dydx)
  {
    real accum = coeff[N];
    dydx = coeff[N] * (real)N;
    for (int i = (int)(N - 1); i > 0; i --) {
      accum = accum * x + coeff[i];
      dydx = dydx * x + coeff[i] * (real)i;
    }

    if (N > 0) {
      accum = accum * x + coeff[0];
    }
    
    return accum;
  }

  bool root(real xmin, real xmax, real epsilon, real threshold, size_t max_iterations, real &x, real &y)
  {
    x = (xmin + xmax)/2.0;
    real yp;
    real delta;

    for (size_t i = 0; i < max_iterations; i ++) {

      y = value_and_dydx(x, yp);

      if (fabs(yp) < epsilon) {
	ERROR("Near zero gradient");
	return false;
      }

      delta = y/yp;
      x -= delta;

      if (x < xmin || x > xmax) {
	ERROR("Root out of range %f %f %f", delta, y, yp);
	return false;
      }
      
      if (fabs(delta) < threshold) {
	// Root found to tolerance
	return true;
      }
    }

    ERROR("Max iterations reached");
    return false;
  }

  Polynomial& operator*=(const real &s)
  {
    for (size_t i = 0; i <= N; i ++) {
      coeff[i] *= s;
    }

    return *this;
  }

  Polynomial& operator/=(const real &s)
  {
    for (size_t i = 0; i <= N; i ++) {
      coeff[i] /= s;
    }

    return *this;
  }

  Polynomial& operator*=(const Polynomial &p)
  {
    Polynomial q(*this);

    size_t new_order = N + p.N;
    resize(new_order + 1);
    for (size_t i = 0; i <= new_order; i ++) {
      coeff[i] = 0.0;
    }
    N = new_order;
    
    for (size_t i = 0; i <= q.N; i ++) {
      for (size_t j = 0; j <= p.N; j ++) {
	coeff[i + j] += q.coeff[i] * p.coeff[j];
      }
    }

    return *this;
  }

  Polynomial& operator+=(const Polynomial &p)
  {
    if (p.N > N) {
      resize(p.N + 1);
      for (size_t i = N + 1; i <= p.N; i ++) {
	coeff[i] = 0.0;
      }
      N = p.N;
    }

    for (size_t i = 0; i <= p.N; i ++) {
      coeff[i] += p.coeff[i];
    }

    return *this;
  }

  Polynomial& operator-=(const Polynomial &p)
  {
    if (p.N > N) {
      resize(p.N + 1);
      for (size_t i = N + 1; i <= p.N; i ++) {
	coeff[i] = 0.0;
      }
      N = p.N;
    }

    for (size_t i = 0; i <= p.N; i ++) {
      coeff[i] -= p.coeff[i];
    }

    return *this;
  }

  friend Polynomial operator*(const Polynomial &p, const real &s)
  {
    Polynomial q(p);

    q *= s;

    return q;
  }

  friend Polynomial operator/(const Polynomial &p, const real &s)
  {
    Polynomial q(p);

    q /= s;

    return q;
  }

  friend Polynomial operator*(const real &s, const Polynomial &p)
  {
    Polynomial q(p);

    q *= s;

    return q;
  }

  friend Polynomial operator*(const Polynomial &p, const Polynomial &q)
  {
    Polynomial r(p);

    r *= q;

    return r;
  }

  friend Polynomial operator+(const Polynomial &p, const Polynomial &q)
  {
    Polynomial r(p);

    r += q;

    return r;
  }

  friend Polynomial operator-(const Polynomial &p, const Polynomial &q)
  {
    Polynomial r(p);

    r -= q;

    return r;
  }

  size_t order() const
  {
    return N;
  }
  
  real operator[](size_t i) const
  {
    if (i <= N) {
      return coeff[i];
    }

    return 0.0;
  }

  Polynomial& operator=(const Polynomial &rhs)
  {
    if (this != &rhs) {
      resize(rhs.size);

      N = rhs.N;
      for (size_t i = 0; i <= N; i ++) {
	coeff[i] = rhs.coeff[i];
      }
    }

    return *this;
  }

  void print(FILE *fp)
  {
    for (size_t i = 0; i <= N; i ++) {

      if (i > 0 && coeff[i] > 0) {
	fprintf(fp, "+");
      }
      
      fprintf(fp, "%.3fx^%d", coeff[i], (int)i);
    }
    fprintf(fp, "\n");
  }
  
private:

  void resize(size_t requiredN)
  {
    size_t new_size = size;

    while (new_size < requiredN) {
      new_size *= 2;
    }

    if (new_size > size) {

      // Need to reallocate

      std::unique_ptr<real[]> new_coeff(new real[new_size]);
      for (size_t i = 0; i < N; i ++) {
	new_coeff[i] = coeff[i];
      }

      coeff = std::move(new_coeff);
      size = new_size;
    }      
  }
  
  // N is order -> no. coeffs = N + 1
  size_t size;
  size_t N;

  // Coeff's are in power order, i.e coeff[0]*x^0 + coeff[1]*x^1 ...
  std::unique_ptr<real[]> coeff;
};

#endif // polynomial_hpp
