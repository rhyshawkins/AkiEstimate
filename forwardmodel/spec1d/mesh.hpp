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
#ifndef mesh_hpp
#define mesh_hpp

#include <vector>

#include "lobattoprojection.hpp"

typedef enum {
  MESHOFFSET_RHO = 0,
  MESHOFFSET_A = 1,
  MESHOFFSET_C = 2,
  MESHOFFSET_F = 3,
  MESHOFFSET_L = 4,
  MESHOFFSET_N = 5
} MeshOffset_t;
  
template
<
  typename real
>
struct MeshParameter
{
  MeshParameter() :
    rho(0.0),
    A(0.0),
    C(0.0),
    F(0.0),
    L(0.0),
    N(0.0)
  {
  }
  
  void zero()
  {
    rho = 0.0;
    A = 0.0;
    C = 0.0;
    F = 0.0;
    L = 0.0;
    N = 0.0;
  }

  void print(FILE *fp)
  {
    fprintf(fp, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
	    (double)rho,
	    (double)A,
	    (double)C,
	    (double)F,
	    (double)L,
	    (double)N);
  }

  MeshParameter& operator+=(const MeshParameter& rhs)
  {
    rho += rhs.rho;
    A += rhs.A;
    C += rhs.C;
    F += rhs.F;
    L += rhs.L;
    N += rhs.N;

    return *this;
  }

  MeshParameter& operator*=(real s)
  {
    rho *= s;
    A *= s;
    C *= s;
    F *= s;
    L *= s;
    N *= s;

    return *this;
  }

  friend MeshParameter operator*(const MeshParameter &rhs, real s)
  {
    MeshParameter p(rhs);
    p *= s;
    return p;
  }
  
  friend MeshParameter operator*(real s, const MeshParameter &rhs)
  {
    MeshParameter p(rhs);
    p *= s;
    return p;
  }

  real rho;
  real A;
  real C;
  real F;
  real L;
  real N;

};

template
<
  typename real
>
struct Amplitude
{
  real U;
  real V;
  real W;

  Amplitude() :
    U(0.0),
    V(0.0),
    W(0.0)
  {
  }

  Amplitude(const Amplitude &a) :
    U(a.U),
    V(a.V),
    W(a.W)
  {
  }

};

template
<
  typename real,
  size_t maxorder
>
struct MeshCell {
  MeshCell(size_t _source_index, real _thickness) :
    source_index(_source_index),
    thickness(_thickness)
  {
  }

  void print(FILE *fp)
  {
    fprintf(fp, "%10.6f\n", (double)thickness);
    for (auto &n : nodes) {
      n.print(fp);
    }
  }

  size_t source_index;
  real thickness;
  size_t order;
  std::array<MeshParameter<real>, maxorder + 1> nodes;

  Spec1DMatrix<real> jacobian;
  real thickness_jacobian;
};

template
<
  typename real,
  size_t maxorder
>
struct CellAmplitude {

  void print(FILE *fp) const
  {
    fprintf(fp, "%d %f\n", (int)order, thickness);
    for (auto &n : nodes) {
      fprintf(fp, "    %10.6f %10.6f %10.6f\n", n.U, n.V, n.W);
    }
  }
  
  real thickness;
  size_t order;
  std::array<Amplitude<real>, maxorder + 1> nodes;
};


template
<
  typename real,
  size_t maxorder
>
class Mesh
{
public:

  Mesh()
  {
    for (size_t i = 0; i <= maxorder; i ++) {
      projection[i] = new LobattoProjection<double, maxorder>(i);
      quadrature[i] = new LobattoQuadrature<double, maxorder>(i);
    }
  }

  ~Mesh()
  {
    for (size_t i = 0; i <= maxorder; i ++) {
      delete projection[i];
      delete quadrature[i];
    }
  }
  
  size_t nodes() const
  {
    size_t c = 1;

    for (auto &cell : cells) {
      c += cell.order;
    }

    return c;
  }

  void print(FILE *fp)
  {
    for (auto &c : cells) {
      fprintf(fp, "Cell\n");
      c.print(fp);
    }
    fprintf(fp, "Basement\n");
    boundary.print(fp);
  }

  MeshParameter<real> interpolate(real z) const
  {
    real dz = z;
    for (auto &c : cells) {

      if (dz < c.thickness) {

	real xi = 2.0*dz/c.thickness - 1.0;

	if (xi < -1.0 || xi > 1.0) {
	  FATAL("xi out of range: %f", xi);
	}
	
	MeshParameter<real> p;

	for (size_t i = 0; i <= c.order; i ++) {

	  p += quadrature[c.order]->cardinal[i].value(xi)*c.nodes[i];

	}

	// if (p.A < 0.0) {
	//   c.print(stderr);
	//   FATAL("-ve A: %f xi %f ", p.A, xi);
	// }
	
	return p;
	
      } else {

	dz -= c.thickness;

      }
    }

    //
    // Return the boundary parameters
    //
    return boundary;
  }
  
  std::array<LobattoProjection<double, maxorder>*, maxorder + 1> projection;
  std::array<LobattoQuadrature<double, maxorder>*, maxorder + 1> quadrature;

  std::vector<MeshCell<real, maxorder>> cells;
  std::vector<size_t> cell_parameter_offsets;
  std::vector<size_t> cell_thickness_reference;

  MeshParameter<real> boundary;
  size_t boundary_parameter_offset;
  Spec1DMatrix<real> boundary_jacobian;
  
};

template
<
  typename real,
  size_t maxorder,
  size_t basementorder = maxorder
>
class MeshAmplitude {
public:
  
  std::vector<CellAmplitude<real, maxorder>> cells;

  CellAmplitude<real, basementorder> boundary;

  void print(FILE *fp)
  {
    for (auto &c: cells) {
      c.print(fp);
    }
    fprintf(fp, "Boundary: \n");
    boundary.print(fp);
  }
   
};

#endif // mesh_hpp

  
