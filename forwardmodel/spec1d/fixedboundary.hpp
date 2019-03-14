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
#ifndef fixedboundary_hpp
#define fixedboundary_hpp

#include "mesh.hpp"

template
<
  typename real,
  size_t maxorder
>
class FixedBoundary {
public:

  FixedBoundary()
  {
  }

  static const char *NAME()
  {
    return "FixedBoundary";
  }

  size_t nparameters() const
  {
    return 0;
  }
  
  void project(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    mesh.boundary.zero();
  }

  void project_gradient(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    mesh.boundary.zero();
    mesh.boundary_parameter_offset = 0;
    mesh.boundary_jacobian.resize(0, 0);
  }
  
  void project(size_t index,
	       real basement,
	       double min_basement,
	       double max_cell_thickness,
	       const MeshParameter<real> &boundary,
	       Mesh<real, maxorder> &mesh) const
  {
    mesh.boundary.zero();

    if (basement < min_basement) {
      //
      // Need to fill in mesh with filler cells
      //
      size_t ncells;
      real cellthickness;
      
      if (max_cell_thickness > 0.0) {

	ncells = ceil((min_basement - basement)/max_cell_thickness);
	cellthickness = (min_basement - basement)/ncells;

      } else {

	ncells = 1;
	cellthickness = (min_basement - basement);

      }

      for (size_t i = 0; i < ncells; i ++) {

	MeshCell<real, maxorder> mcell(index + i, cellthickness);
	mcell.order = maxorder;

	for (size_t k = 0; k <= maxorder; k ++) {

	  mcell.nodes[k].rho = boundary.rho;
	  mcell.nodes[k].A = boundary.A;
	  mcell.nodes[k].C = boundary.C;
	  mcell.nodes[k].F = boundary.F;
	  mcell.nodes[k].L = boundary.L;
	  mcell.nodes[k].N = boundary.N;
	    
	}

	mesh.cells.push_back(mcell);
      }
    }
  }

  void print(FILE *fp) const
  {
    fprintf(fp, "FixedBoundary\n");
  }

  bool read(FILE *fp)
  {
    int order;
    double thickness;
    
    if (fscanf(fp, "%d %lf\n", &order, &thickness) != 2) {
      ERROR("Failed to read cell parameters");
      return false;
    }

    if (order != 0 || thickness != 0.0) {
      ERROR("Expected terminating cell parameters, got %d %f", order, thickness);
      return false;
    }
      
    return true;
  }

  bool save(FILE *fp) const
  {
    return true;
  }

  int encode_size()
  {
    return 0;
  }
  
  int encode(char *buffer, int &buffer_offset, int buffer_size)
  {
    return 0;
  }
  
  int decode(const char *buffer, int &buffer_offset, int buffer_size)
  {
    return 0;
  }

  static bool is_fixed()
  {
    return true;
  }

  std::array<real, 0> parameters;
};

#endif // fixedboundary_hpp
