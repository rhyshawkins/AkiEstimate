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
#ifndef anisotropicrhovsxivpvshalfspace_hpp
#define anisotropicrhovsxivpvshalfspace_hpp

#include "mesh.hpp"

#include "halfspace.hpp"

template
<
  typename real,
  size_t maxorder
>
class AnisotropicRhoVsXiVpVsHalfspace : public Halfspace<real, 4> {
public:

  enum {
    pRHO = 0,
    pVS,
    pXI,
    pVPVS
  };
  
  AnisotropicRhoVsXiVpVsHalfspace()
  {
  }

  AnisotropicRhoVsXiVpVsHalfspace(real rho,
				  real Vs,
				  real Xi,
				  real VpVs) 
  {
    (*this)[pRHO] = rho;
    (*this)[pVS] = Vs;
    (*this)[pXI] = Xi;
    (*this)[pVPVS] = VpVs;
  }

  static const char *NAME()
  {
    return "AnisotropicRhoVsXiVpVsHalfspace";
  }

  void project(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    real vp = (*this)[pVS] * (*this)[pVPVS];
    mesh.boundary.rho = (*this)[pRHO];
    mesh.boundary.A = (*this)[pRHO] * vp*vp;
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = 3.0*(*this)[pRHO] * (*this)[pVS] * (*this)[pVS]/(2.0 + (*this)[pXI]);
    mesh.boundary.N = 3.0*(*this)[pRHO] * (*this)[pVS] * (*this)[pVS]/(2.0/(*this)[pXI] + 1.0);
    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;
  }

  void project_gradient(size_t index, real depth, Mesh<real, maxorder> &mesh) const
  {
    //
    // rho, vs, xi, vpvs
    //
    real rho = (*this)[pRHO];
    real vs = (*this)[pVS];
    real xi = (*this)[pXI];
    real vpvs = (*this)[pVPVS];
    real vp = vs * vpvs;
    
    real vp2 = vp * vp;
    real vs2 = vs * vs;

    mesh.boundary.rho = rho;
    mesh.boundary.A = (*this)[pRHO] * vp * vp;
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = 3.0 * rho * vs2/(2.0 + xi);
    mesh.boundary.N = 3.0 * rho * vs2/(2.0/xi + 1.0);

    mesh.boundary.F = mesh.boundary.A - 2.0*mesh.boundary.L;

    mesh.boundary_parameter_offset = mesh.cell_parameter_offsets[mesh.cell_parameter_offsets.size() - 1];
    mesh.boundary_jacobian.resize(4, 6);
    mesh.boundary_jacobian.setZero();

    // d/drho
    mesh.boundary_jacobian(0, 0) = 1.0;
    mesh.boundary_jacobian(0, 1) = vp2;
    mesh.boundary_jacobian(0, 2) = vp2;
    mesh.boundary_jacobian(0, 3) = vp2 - 6.0*vs2/(2.0 + xi);
    mesh.boundary_jacobian(0, 4) = 3.0*vs2/(2.0 + xi);
    mesh.boundary_jacobian(0, 5) = 3.0*vs2/(2.0/xi + 1.0);

    // d/dvs
    mesh.boundary_jacobian(1, 0) = 0.0;
    mesh.boundary_jacobian(1, 1) = 2.0*rho*vs*vpvs*vpvs;
    mesh.boundary_jacobian(1, 2) =
      mesh.boundary_jacobian(1, 1) -
      2.0*mesh.boundary_jacobian(1, 4);
    mesh.boundary_jacobian(1, 4) = 6.0*rho*vs/(2.0 + xi);
    mesh.boundary_jacobian(1, 5) = 6.0*rho*vs/(2.0/xi + 1.0);

    mesh.boundary_jacobian(1, 3) = 
      mesh.boundary_jacobian(1, 1) -
      2.0*mesh.boundary_jacobian(1, 4);
    
    // d/dxi
    mesh.boundary_jacobian(2, 0) = 0.0;
    mesh.boundary_jacobian(2, 1) = 0.0;
    mesh.boundary_jacobian(2, 2) = 0.0;
    mesh.boundary_jacobian(2, 4) = -3.0 * rho * vs2 / ((2.0 + xi) * (2.0 + xi));
    mesh.boundary_jacobian(2, 5) = 6.0 * rho * vs2 / ((2.0 + xi) * (2.0 + xi));
    mesh.boundary_jacobian(2, 3) = -2.0*mesh.boundary_jacobian(2, 4);

    // d/dvpvs
    mesh.boundary_jacobian(3, 0) = 0.0;
    mesh.boundary_jacobian(3, 1) = 2.0 * vp2 * vpvs;
    mesh.boundary_jacobian(3, 2) = 2.0 * vp2 * vpvs;
    mesh.boundary_jacobian(3, 3) = mesh.boundary_jacobian(3, 1);
    mesh.boundary_jacobian(3, 4) = 0.0;
    mesh.boundary_jacobian(3, 5) = 0.0;
  }

  void project(size_t index,
	       real basement,
	       double min_basement,
	       double max_cell_thickness,
	       const MeshParameter<real> &boundary,
	       Mesh<real, maxorder> &mesh) const
  {
    real vp = (*this)[pVS] * (*this)[pVPVS];
    mesh.boundary.rho = (*this)[pRHO];
    mesh.boundary.A = (*this)[pRHO] * vp*vp;
    mesh.boundary.C = mesh.boundary.A;
    mesh.boundary.L = 3.0*(*this)[pRHO] * (*this)[pVS] * (*this)[pVS]/(2.0 + (*this)[pXI]);
    mesh.boundary.N = 3.0*(*this)[pRHO] * (*this)[pVS] * (*this)[pVS]/(2.0/(*this)[pXI] + 1.0);
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

#endif // anisotropicrhovsxivpvshalfspace_hpp
