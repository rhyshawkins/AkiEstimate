#pragma once
#ifndef spec1dmatrix_hpp
#define spec1dmatrix_hpp

#include <string.h>

#include "logging.hpp"

#define SPEC1DMATRIX_BOUNDSCHECK

template
<
  typename real
>
class Spec1DMatrix {
public:

  Spec1DMatrix() :
    r(0),
    c(0),
    s(0),
    d(nullptr)
  {
  }

  Spec1DMatrix(const Spec1DMatrix &rhs) :
    r(rhs.r),
    c(rhs.c),
    s(rhs.s),
    d(nullptr)
  {
    if (s > 0) {
      d = new real[s];
      for (int i = 0; i < s; i ++) {
	d[i] = rhs.d[i];
      }
    }
  }
  
  ~Spec1DMatrix()
  {
    delete [] d;
  }

  int rows() const
  {
    return r;
  }

  int cols() const
  {
    return c;
  }

  real *data()
  {
    return d;
  }

  const real *data() const
  {
    return d;
  }

  void resize(int _rows, int _cols)
  {
    r = _rows;
    c = _cols;

    int newsize = r * c;
    
    checksize(newsize);
  }
  
  void setZero()
  {
    if (s > 0) {
      for (int i = 0; i < s; i ++) {
	d[i] = 0.0;
      }
      // memset(d, 0, sizeof(real) * s);
    }
  }

  void setIdentity(int N)
  {
    resize(N, N);
    setZero();
    for (int i = 0; i < N; i ++) {
      d[i * N + i] = 1.0;
    }
  }

  real norm() const
  {
    int t = r * c;
    real sum = 0.0;
    for (int i = 0; i < t; i ++) {
      sum += d[i]*d[i];
    }

    return sqrt(sum);
  }

  
  real &operator()(int i, int j)
  {
    // Matrix stored in column major order (fortran style)
#ifdef SPEC1DMATRIX_BOUNDSCHECK
    if (i < 0 || i >= r) {
      FATAL("Row out of range: %d (%d)", i, r);
    }
    if (j < 0 || j >= c) {
      FATAL("Col out of range: %d (%d)", j, c);
    }
#endif // SPEC1DMATRIX_BOUNDSCHECK
    return d[j * r + i];
  }
  
  const real &operator()(int i, int j) const
  {
#ifdef SPEC1DMATRIX_BOUNDSCHECK
    if (i < 0 || i >= r) {
      FATAL("Row out of range: %d (%d)", i, r);
    }
    if (j < 0 || j >= c) {
      FATAL("Col out of range: %d (%d)", j, c);
    }
#endif // SPEC1DMATRIX_BOUNDSCHECK
    return d[j * r + i];
  }

  Spec1DMatrix &operator=(const Spec1DMatrix &rhs)
  {
    if (&rhs != this) {
      resize(rhs.r, rhs.c);

      int n = rhs.r * rhs.c;
      for (int i = 0; i < n; i ++) {
	d[i] = rhs.d[i];
      }
    }

    return *this;
  }
  
  real *col(int j)
  {
#ifdef SPEC1DMATRIX_BOUNDSCHECK
    if (j < 0 || j >= c) {
      FATAL("Col out of range: %d (%d)", j, c);
    }
#endif // SPEC1DMATRIX_BOUNDSCHECK
    return d + j * r;
  }
  
  const real *col(int j) const
  {
#ifdef SPEC1DMATRIX_BOUNDSCHECK
    if (j < 0 || j >= c) {
      FATAL("Col out of range: %d (%d)", j, c);
    }
#endif // SPEC1DMATRIX_BOUNDSCHECK
    return d + j * r;
  }


  Spec1DMatrix &operator+=(const Spec1DMatrix &rhs)
  {
    if (rhs.r != r || rhs.c != c) {
      FATAL("Array mismatch");
    }
    
    int t = r*c;
    for (int i = 0; i < t; i ++) {
      d[i] += rhs.d[i];
    }

    return *this;
  }
  
  Spec1DMatrix &operator-=(const Spec1DMatrix &rhs)
  {
    if (rhs.r != r || rhs.c != c) {
      FATAL("Array mismatch");
    }
    
    int t = r*c;
    for (int i = 0; i < t; i ++) {
      d[i] -= rhs.d[i];
    }

    return *this;
  }

  friend Spec1DMatrix operator*(const Spec1DMatrix &a, real s)
  {
    Spec1DMatrix r(a);
    r *= s;
    return r;
  }
  
  friend Spec1DMatrix operator*(real s, const Spec1DMatrix &a)
  {
    Spec1DMatrix r(a);
    r *= s;
    return r;
  }

  Spec1DMatrix &operator*=(real s)
  {
    int t = r*c;
    for (int i = 0; i < t; i ++) {
      d[i] *= s;
    }

    return *this;
  }
  
  Spec1DMatrix &operator/=(real s)
  {
    int t = r*c;
    for (int i = 0; i < t; i ++) {
      d[i] /= s;
    }

    return *this;
  }

private:

  void checksize(int required_size)
  {
    if (required_size > s) {
      while (required_size > s) {
	if (s == 0) {
	  s = 1024;
	} else {
	  s = 2*s;
	}
      }

      if (d != nullptr) {
	delete [] d;
      }
      
      d = new real[s];
    }
  }

  int r;
  int c;
  int s;

  real *d;

};


#endif // spec1dmatrix_hpp
  

    
