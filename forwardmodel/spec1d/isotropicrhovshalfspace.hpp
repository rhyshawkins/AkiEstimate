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
#ifndef isotropicrhovshalfspace_hpp
#define isotropicrhovshalfspace_hpp

#include "mesh.hpp"

template
<
  typename real,
  size_t maxorder,
  size_t vpvs_numerator = 17,
  size_t vpvs_denominator = 10
>
class IsotropicRhoVsHalfspace {
public:

  IsotropicRhoVsHalfspace()
  {
  }

  IsotropicRhoVsHalfspace(real rho,
			  real Vs)
  {
    parameters[0] = rho;
    parameters[1] = Vs;
  }

  void project(real depth, Mesh<real, maxorder> &mesh) const
  {
    real vpvs = (real)vpvs_numerator/(real)vpvs_denominator;

    mesh.boundary.rho = parameters[0];
    mesh.boundary.A = parameters[0] * (vpvs * parameters[1] * vpvs * parameters[1]);
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = parameters[0] * parameters[1] * parameters[1];
    mesh.boundary.N = mesh.boundary.L;
    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;
  }

  void project(real basement,
	       double min_basement,
	       double max_cell_thickness,
	       const MeshParameter<real> &boundary,
	       Mesh<real, maxorder> &mesh) const
  {
    real vpvs = (real)vpvs_numerator/(real)vpvs_denominator;

    mesh.boundary.rho = parameters[0];
    mesh.boundary.A = parameters[0] * (vpvs * parameters[1] * vpvs * parameters[1]);
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = parameters[0] * parameters[1] * parameters[1];
    mesh.boundary.N = mesh.boundary.L;
    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;

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

	MeshCell<real, maxorder> mcell(cellthickness);

	for (size_t k = 0; k <= maxorder; k ++) {

	  mcell.nodes[k].rho = mesh.boundary.rho;
	  mcell.nodes[k].A = mesh.boundary.A;
	  mcell.nodes[k].C = mesh.boundary.C;
	  mcell.nodes[k].F = mesh.boundary.F;
	  mcell.nodes[k].L = mesh.boundary.L;
	  mcell.nodes[k].N = mesh.boundary.N;

	}

	mesh.cells.push_back(mcell);
      }
    }
  }

  bool read(FILE *fp)
  {
    int nparameters;
    double thickness;
    double rho, vs;
    
    if (fscanf(fp, "%d %lf\n", &nparameters, &thickness) != 2) {
      ERROR("Failed to read cell parameters");
      return false;
    }

    if (nparameters != 2 || thickness != 0.0) {
      ERROR("Expected terminating cell parameters, got %d %f", nparameters, thickness);
      return false;
    }
      
    if (fscanf(fp, "%lf %lf\n", &rho, &vs) != 2) {
      ERROR("Failed to read parameters");
      return false;
    }

    parameters[0] = rho;
    parameters[1] = vs;

    return true;
  }
  
  bool save(FILE *fp) const
  {
    fprintf(fp, "%15.9f %15.9f\n", (double)parameters[0], (double)parameters[1]);
    return true;
  }

  void print(FILE *fp) const
  {
    fprintf(fp, "IsotropicRhoVsHalfspace\n%15.9f %15.9f\n", (double)parameters[0], (double)parameters[1]);
  }

  int encode_size()
  {
    return 2*sizeof(real);
  }
  
  int encode(char *buffer, int &buffer_offset, int buffer_size)
  {
    if (::encode<real>(parameters[0], buffer, buffer_offset, buffer_size) < 0) {
      return -1;
    }
    if (::encode<real>(parameters[1], buffer, buffer_offset, buffer_size) < 0) {
      return -1;
    }

    return 2*sizeof(real);
  }

  int decode(const char *buffer, int &buffer_offset, int buffer_size)
  {
    if (::decode<real>(parameters[0], buffer, buffer_offset, buffer_size) < 0) {
      return -1;
    }
    if (::decode<real>(parameters[1], buffer, buffer_offset, buffer_size) < 0) {
      return -1;
    }

    return 2*sizeof(real);
  }

  static bool is_fixed()
  {
    return false;
  }

  // Order is Rho, Vs
  std::array<real, 2> parameters;
  
};

#endif // isotropicrhovshalfspace_hpp
