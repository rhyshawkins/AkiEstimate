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
#ifndef isotropicvshalfspace_hpp
#define isotropicvshalfspace_hpp

#include "mesh.hpp"

#include "halfspace.hpp"

template
<
  typename real,
  typename empiricalmodel,
  size_t maxorder
>
class IsotropicVsHalfspace : public Halfspace<real, 1> {
public:

  IsotropicVsHalfspace() :
    Halfspace<real, 1>()
  {
  }

  IsotropicVsHalfspace(real Vs)
  {
    (*this)[0] = Vs;
  }

  static const char *NAME()
  {
    return "IsotropicVsHalfspace";
  }

  void project(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    real vp, rho;
    
    empiricalmodel::compute((*this)[0], rho, vp);
    
    mesh.boundary.rho = rho;
    mesh.boundary.A = rho * (vp * vp);
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = rho * (*this)[0] * (*this)[0];
    mesh.boundary.N = mesh.boundary.L;
    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;
  }

  void project_gradient(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    real vp, rho, dvpdvs, drhodvs;
    
    empiricalmodel::compute_gradient((*this)[0], rho, vp, drhodvs, dvpdvs);
    
    mesh.boundary.rho = rho;
    mesh.boundary.A = rho * (vp * vp);
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = rho * (*this)[0] * (*this)[0];
    mesh.boundary.N = mesh.boundary.L;
    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;

    mesh.boundary_parameter_offset = mesh.cell_parameter_offsets[mesh.cell_parameter_offsets.size() - 1];
    mesh.boundary_jacobian.resize(1, 6);
    mesh.boundary_jacobian.setZero();

    mesh.boundary_jacobian(0, 0) = drhodvs;
    mesh.boundary_jacobian(0, 1) = drhodvs * vp * vp + rho * 2.0 * vp * dvpdvs;
    mesh.boundary_jacobian(0, 2) = mesh.boundary_jacobian(0, 1);
    mesh.boundary_jacobian(0, 3) = drhodvs*vp*vp + rho*2.0*vp*dvpdvs -
      2.0*(drhodvs*(*this)[0] + rho*2.0)*(*this)[0];
    mesh.boundary_jacobian(0, 4) = (drhodvs*(*this)[0] + rho*2.0)*(*this)[0];
    mesh.boundary_jacobian(0, 5) = mesh.boundary_jacobian(0, 4);
  }

  void project(size_t index, real basement,
	       double min_basement,
	       double max_cell_thickness,
	       const MeshParameter<real> &boundary,
	       Mesh<real, maxorder> &mesh) const
  {
    real vp, rho;

    empiricalmodel::compute((*this)[0], rho, vp);
    
    mesh.boundary.rho = rho;
    mesh.boundary.A = rho * (vp * vp);
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = rho * (*this)[0] * (*this)[0];
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


};

#endif // isotropicvshalfspace_hpp
