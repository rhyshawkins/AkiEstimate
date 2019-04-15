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
#ifndef isotropicrhovpvshalfspace_hpp
#define isotropicrhovpvshalfspace_hpp

#include "mesh.hpp"

template
<
  typename real,
  size_t maxorder
>
class IsotropicRhoVpVsHalfspace {
public:

  IsotropicRhoVpVsHalfspace()
  {
  }

  IsotropicRhoVpVsHalfspace(real rho,
			    real Vp,
			    real Vs)
  {
    parameters[0] = rho;
    parameters[1] = Vp;
    parameters[2] = Vs;
  }

  size_t nparameters() const
  {
    return 0;
  }

  void project(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    mesh.boundary.rho = parameters[0];
    mesh.boundary.A = parameters[0] * parameters[1] * parameters[1];
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = parameters[0] * parameters[2] * parameters[2];
    mesh.boundary.N = mesh.boundary.L;
    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;
  }

  void project_gradient(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    real vp2 = parameters[1]*parameters[1];
    real vs2 = parameters[2]*parameters[2];

    mesh.boundary.rho = parameters[0];
    mesh.boundary.A = parameters[0] * vp2;
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = parameters[0] * vs2;
    mesh.boundary.N = mesh.boundary.L;

    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;

    mesh.boundary_parameter_offset = mesh.cell_parameter_offsets[mesh.cell_parameter_offsets.size() - 1];
    mesh.boundary_jacobian.resize(3, 6);
    mesh.boundary_jacobian.setZero();

    // d/drho
    mesh.boundary_jacobian(0, 0) = 1.0;
    mesh.boundary_jacobian(0, 1) = vp2;
    mesh.boundary_jacobian(0, 2) = vp2;
    mesh.boundary_jacobian(0, 3) = vp2 - 2.0*vs2;
    mesh.boundary_jacobian(0, 4) = vs2;
    mesh.boundary_jacobian(0, 5) = vs2;

    // d/dvp
    mesh.boundary_jacobian(1, 0) = 0.0;
    mesh.boundary_jacobian(1, 1) = 2.0*parameters[0]*parameters[1];
    mesh.boundary_jacobian(1, 2) = mesh.boundary_jacobian(1, 1);
    mesh.boundary_jacobian(1, 3) = mesh.boundary_jacobian(1, 1);
    mesh.boundary_jacobian(1, 4) = 0.0;
    mesh.boundary_jacobian(1, 5) = 0.0;
    
    // d/dvs
    mesh.boundary_jacobian(2, 0) = 0.0;
    mesh.boundary_jacobian(2, 1) = 0.0;
    mesh.boundary_jacobian(2, 2) = 0.0;
    mesh.boundary_jacobian(2, 3) = -4.0 * parameters[0] * parameters[2];
    mesh.boundary_jacobian(2, 4) = 2.0 * parameters[0] * parameters[2];
    mesh.boundary_jacobian(2, 5) = mesh.boundary_jacobian(2, 4);
  }

  void project(size_t index,
	       real basement,
	       double min_basement,
	       double max_cell_thickness,
	       const MeshParameter<real> &boundary,
	       Mesh<real, maxorder> &mesh) const
  {
    mesh.boundary.rho = parameters[0];
    mesh.boundary.A = parameters[0] * parameters[1] * parameters[1];
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = parameters[0] * parameters[2] * parameters[2];
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

	MeshCell<real, maxorder> mcell(index + i, cellthickness);
	mcell.order = maxorder;
	
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
    
    if (fscanf(fp, "%d %lf\n", &nparameters, &thickness) != 2) {
      ERROR("Failed to read cell parameters");
      return false;
    }

    if (nparameters != 3 || thickness != 0.0) {
      ERROR("Expected terminating cell parameters, got %d %f", nparameters, thickness);
      return false;
    }

    int order;
    double rho, vp, vs;

    if (fscanf(fp, "%d %lf\n", &order, &rho) != 2) {
      ERROR("Failed to read density");
      return false;
    }

    if (order != 0) {
      ERROR("Invalid order for density %d", order);
      return false;
    }

    if (fscanf(fp, "%d %lf\n", &order, &vp) != 2) {
      ERROR("Failed to read Vp");
      return false;
    }

    if (order != 0) {
      ERROR("Invalid order for Vs %d", order);
      return false;
    }
    
    if (fscanf(fp, "%d %lf\n", &order, &vs) != 2) {
      ERROR("Failed to read Vs");
      return false;
    }

    if (order != 0) {
      ERROR("Invalid order for Vs %d", order);
      return false;
    }

    parameters[0] = rho;
    parameters[1] = vp;
    parameters[2] = vs;

    return true;
  }

  static bool is_fixed()
  {
    return false;
  }

  // Order is Rho, Vp, Vs
  std::array<real, 3> parameters;
  
};

#endif // isotropicrhovpvshalfspace_hpp
