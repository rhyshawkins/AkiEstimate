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
#ifndef lovematrices_hpp
#define lovematrices_hpp

#include <vector>

#include "lobattoquadrature.hpp"
#include "laguerrequadrature.hpp"

#include "mesh.hpp"

#ifdef USE_DGGEVPP
#include "dggev.hpp"
#else 
#include "generalisedeigenproblem.hpp"
#endif // USE_DGGEVPP

#include "specializedeigenproblem.hpp"
#include "spec1dmatrix.hpp"

template
<
  typename real,
  size_t maxorder,
  size_t maxboundaryorder = maxorder 
>
class LoveMatrices {
public:

  LoveMatrices() :
    size(0)
  {
    for (size_t i = 0; i <= maxorder; i ++) {
      Lobatto[i] = new LobattoQuadrature<double, maxorder>(i);
    }
    
    for (size_t i = 0; i <= maxboundaryorder; i ++) {
      Laguerre[i] = new LaguerreQuadrature<double, maxboundaryorder>(i);
    }
  }

  ~LoveMatrices()
  {
    for (size_t i = 0; i <= maxorder; i ++) {
      delete Lobatto[i];
    }
    
    for (size_t i = 0; i <= maxboundaryorder; i ++) {
      delete Laguerre[i];
    }
  }
  
  void recompute(const Mesh<real, maxorder> &mesh, size_t boundaryorder, real scale = 1.0)
  {
    if (boundaryorder > maxboundaryorder) {
      FATAL("Boundary order out of range: %d > %d", (int)boundaryorder, (int)maxboundaryorder);
    }
    
    size_t cell_nodes = mesh.nodes();

    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre boundary (extra order - 1 nodes)
      //
      size = cell_nodes + boundaryorder;

    } else {
      //
      // Fixed boundary (last node is zero and removed from system of equations)
      //
      size = cell_nodes - 1;

    }

    laguerrescale = scale;
    
    A.resize(size, size);
    B.resize(size, size);
    Bs.resize(size, size);
    
    C.resize(size, size);
    D.resize(size, size);

    dB.resize(size, size);
    dD.resize(size, size);

    computeA(mesh, boundaryorder);
    computeB(mesh, boundaryorder);
    computeC(mesh, boundaryorder);

    l2values.resize(mesh.cells.size() + 1, 2);
  }

  size_t postcomputegradient(const Mesh<real, maxorder> &mesh,
			     size_t boundaryorder,
			     Spec1DMatrix<real> &v,
			     size_t &nbasecells)
  {
    size_t nparameters = 0;
    for (auto &c: mesh.cells) {
      nparameters += c.jacobian.rows();
    }
    if (nparameters == 0) {
      FATAL("Jacobians not computed");
    }
    
    if (mesh.boundary.rho > 0.0) {
      nparameters += mesh.boundary_jacobian.rows();
    }

    if (mesh.cell_thickness_reference.size() == 0) {
      FATAL("Thickness references not set");
    }

    nbasecells = mesh.cell_thickness_reference[mesh.cell_thickness_reference.size() - 1] + 1;

    computeAGradient(nparameters, nbasecells, mesh, boundaryorder, v);
    computeBGradient(nparameters, nbasecells, mesh, boundaryorder, v);
    computeCGradient(nparameters, nbasecells, mesh, boundaryorder, v);

    return nparameters;
  }

  real solve_fundamental_gep(real omega, real &normA, real &normB, real &normC)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i));
	Bs(j, i) = B(j, i);
      }
    }

    if (!GEP(D, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalised eigen problem");
      return 0.0;
    }
    
    real fundamental = 0.0;

    u.resize(size, 1);
    v.resize(size, 1);

    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0 && k2 > fundamental) {
	  fundamental = k2;

	  for (size_t j = 0; j < size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    if (fundamental > 0.0) {
      
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      for (size_t j = 0; j < size; j ++) {
	normA += u(j, 0) * A(j, j) * v(j, 0);
	normB += u(j, 0) * B(j, j) * v(j, 0);

	real c = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  c += C(j, i) * v(i, 0);
	}
	normC += u(j, 0) * c;
      }

      

    }
    return sqrt(fundamental);
  }

  real solve_fundamental_gep_continuity(const Mesh<real, maxorder> &mesh,
					size_t boundaryorder,
					real omega,
					real &normA, real &normB, real &normC,
					real &normL2)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i));
	Bs(j, i) = B(j, i);
      }
    }

    if (!GEP(D, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalised eigen problem");
      return 0.0;
    }
    
    real fundamental = 0.0;

    u.resize(size, 1);
    v.resize(size, 1);

    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0 && k2 > fundamental) {
	  fundamental = k2;

	  for (size_t j = 0; j < size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    if (fundamental > 0.0) {
      
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      for (size_t j = 0; j < size; j ++) {
	normA += u(j, 0) * A(j, j) * v(j, 0);
	normB += u(j, 0) * B(j, j) * v(j, 0);

	real c = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  c += C(j, i) * v(i, 0);
	}
	normC += u(j, 0) * c;
      }

      size_t poffset = 0;
      size_t ci = 0;
      for (auto &c : mesh.cells) {

	real dV0 = 0.0;
	real dVn = 0.0;
	
	for (size_t i = 0; i <= c.order; i ++) {

	  dV0 += Lobatto[c.order]->cardinal[i].dydz(-1.0) * v(poffset + i, 0) * 2.0/c.thickness;
	  dVn += Lobatto[c.order]->cardinal[i].dydz(1.0) * v(poffset + i, 0) * 2.0/c.thickness;

	}

	l2values(ci, 0) = dV0;
	l2values(ci, 1) = dVn;
	
	poffset += (c.order);
	ci ++;
	
      }

      if (mesh.boundary.rho > 0.0) {
	real dV0;
	for (size_t i = 0; i <= boundaryorder; i ++) {
	  dV0 += Laguerre[boundaryorder]->cardinal[i].dydz(-1.0) * v(poffset + i, 0) * laguerrescale;
	}

	l2values(ci, 0) = dV0;
	l2values(ci, 1) = 0.0;
      }

      normL2 = 0.0;
      for (size_t i = 0; i < ci; i ++) {
	printf("%2d %16.9e %16.9e : %16.9e\n", (int)i, l2values(i, 0), l2values(i, 1));
      }
      
    }
    return sqrt(fundamental);
  }
  
  real solve_fundamental_sep(real omega, real &normA, real &normB, real &normC)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i));
	Bs(j, i) = B(j, i);
      }
    }

    if (!SpecializedEigenProblem(D, Bs, work, eu, ev, lambda, Q, Z)) {
      ERROR("Failed to compute generalised eigen problem");
      return 0.0;
    }
    
    real fundamental = 0.0;

    u.resize(size, 1);
    v.resize(size, 1);

    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0 && k2 > fundamental) {
	  fundamental = k2;

	  for (size_t j = 0; j < size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    if (fundamental > 0.0) {
      
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      for (size_t j = 0; j < size; j ++) {
	normA += u(j, 0) * A(j, j) * v(j, 0);
	normB += u(j, 0) * B(j, j) * v(j, 0);

	real c = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  c += C(j, i) * v(i, 0);
	}
	normC += u(j, 0) * c;
      }

    }
    return sqrt(fundamental);
  }

  real solve_fundamental_gradient(const Mesh<real, maxorder> &mesh,
				  size_t boundaryorder,
				  real omega,
				  Spec1DMatrix<real> &dkdp,
				  Spec1DMatrix<real> &dUdp,
				  real &normA,
				  real &normB,
				  real &normC)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i));
	Bs(j, i) = B(j, i);
      }
    }

    if (!GEP(D, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalised eigen problem");
      return 0.0;
    }
    
    real fundamental = 0.0;
    v.resize(size, 1);
    u.resize(size, 1);

    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0 && k2 > fundamental) {
	  fundamental = k2;

	  for (size_t j = 0; j < size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    double unorm = 0.0;
    double vnorm = 0.0;
    for (size_t j = 0; j < size; j ++) {
      unorm += u(j, 0)*u(j, 0);
      vnorm += v(j, 0)*v(j, 0);
    }
    unorm = sqrt(unorm);
    vnorm = sqrt(vnorm);

    if (v(0, 0) > 0.0) {
      for (size_t j = 0; j < size; j ++) {
	u(j, 0) /= unorm;
	v(j, 0) /= vnorm;
      }
    } else {
      for (size_t j = 0; j < size; j ++) {
	u(j, 0) /= -unorm;
	v(j, 0) /= -vnorm;
      }
    }      
    
    real k = 0.0;
    
    if (fundamental > 0.0) {
      k = sqrt(fundamental);

      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      for (size_t j = 0; j < size; j ++) {
	normA += v(j, 0) * A(j, j) * v(j, 0);
	normB += v(j, 0) * B(j, j) * v(j, 0);
	real c = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  c += C(j, i) * v(i, 0);
	}
	normC += v(j, 0) * c;
      }

      size_t nbasecells;
      size_t nparameters = postcomputegradient(mesh,
					       boundaryorder,
					       v,
					       nbasecells);

      dkdp.resize(nparameters, 1);
      dkdp.setZero();

      dUdp.resize(nparameters, 1);
      dUdp.setZero();

      // FILE *fp = fopen("loveU.txt", "a");
      
      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;
	
	for (size_t i = 0; i < size; i ++) {
	  udAv += v(i, 0) * dAv(i, j);
	  udBv += v(i, 0) * dBv(i, j);
	  udCv += v(i, 0) * dCv(i, j);
	}

	dkdp(j, 0) = ((o2*udAv - udCv) - fundamental*udBv)/(2.0 * k * normB);

	dUdp(j, 0) = (normA*(dkdp(j, 0)*normB + k*udBv) - k*normB*udAv)/(omega*normA*normA);
	//dUdp(j, 0) = normB/(omega * normA) * dkdp(j, 0);

	// dUdp(j, 0) = dkdp(j, 0)*normB + k*udBv; // For testing d/dp Bk
	// dUdp(j, 0) = omega*udAv; // For testing d/dp omega A
	// dUdp(j, 0) = -udAv/(omega*normA*normA); // For testing d/dp 1/(omega A)
	
	if (j == 0) {
	   // fprintf(fp, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	   // 	   normA, normB, normC, dkdp(j, 0), k, udAv, udBv, udCv, omega);

	  // printf("%10.3e %10.3e %10.3e %10.3e\n",
	  // 	 normA*(dkdp(j, 0)*normB),
	  // 	 normA*k*udBv,
	  // 	 k*normB*udAv,
	  // 	 dkdp(j, 0));
	}
      }

      // fclose(fp);
    }

    return k;
  }

  real solve_fundamental_gradient_thickness(const Mesh<real, maxorder> &mesh,
					    size_t boundaryorder,
					    real omega,
					    Spec1DMatrix<real> &dkdp,
					    Spec1DMatrix<real> &dUdp,
					    Spec1DMatrix<real> &dkdT,
					    real &normA,
					    real &normB,
					    real &normC)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i));
	Bs(j, i) = B(j, i);
      }
    }

    if (!GEP(D, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalised eigen problem");
      return 0.0;
    }
    
    real fundamental = 0.0;
    v.resize(size, 1);
    u.resize(size, 1);

    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0 && k2 > fundamental) {
	  fundamental = k2;

	  for (size_t j = 0; j < size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    double unorm = 0.0;
    double vnorm = 0.0;
    for (size_t j = 0; j < size; j ++) {
      unorm += u(j, 0)*u(j, 0);
      vnorm += v(j, 0)*v(j, 0);
    }
    unorm = sqrt(unorm);
    vnorm = sqrt(vnorm);

    if (v(0, 0) > 0.0) {
      for (size_t j = 0; j < size; j ++) {
	u(j, 0) /= unorm;
	v(j, 0) /= vnorm;
      }
    } else {
      for (size_t j = 0; j < size; j ++) {
	u(j, 0) /= -unorm;
	v(j, 0) /= -vnorm;
      }
    }      
    
    real k = 0.0;
    
    if (fundamental > 0.0) {
      k = sqrt(fundamental);

      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      for (size_t j = 0; j < size; j ++) {
	normA += v(j, 0) * A(j, j) * v(j, 0);
	normB += v(j, 0) * B(j, j) * v(j, 0);
	real c = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  c += C(j, i) * v(i, 0);
	}
	normC += v(j, 0) * c;
      }

      size_t nbasecells;
      size_t nparameters = postcomputegradient(mesh,
					       boundaryorder,
					       v,
					       nbasecells);

      dkdp.resize(nparameters, 1);
      dkdp.setZero();

      dUdp.resize(nparameters, 1);
      dUdp.setZero();

      dkdT.resize(nbasecells, 1);
      dkdT.setZero();

      // FILE *fp = fopen("loveU.txt", "a");
      
      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;

	for (size_t i = 0; i < size; i ++) {
	  udAv += v(i, 0) * dAv(i, j);
	  udBv += v(i, 0) * dBv(i, j);
	  udCv += v(i, 0) * dCv(i, j);
	}

	dkdp(j, 0) = ((o2*udAv - udCv) - fundamental*udBv)/(2.0 * k * normB);

	dUdp(j, 0) = (normA*(dkdp(j, 0)*normB + k*udBv) - k*normB*udAv)/(omega*normA*normA);
	//dUdp(j, 0) = normB/(omega * normA) * dkdp(j, 0);

	// dUdp(j, 0) = dkdp(j, 0)*normB + k*udBv; // For testing d/dp Bk
	// dUdp(j, 0) = omega*udAv; // For testing d/dp omega A
	// dUdp(j, 0) = -udAv/(omega*normA*normA); // For testing d/dp 1/(omega A)
	
	if (j == 0) {
	   // fprintf(fp, "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
	   // 	   normA, normB, normC, dkdp(j, 0), k, udAv, udBv, udCv, omega);

	  // printf("%10.3e %10.3e %10.3e %10.3e\n",
	  // 	 normA*(dkdp(j, 0)*normB),
	  // 	 normA*k*udBv,
	  // 	 k*normB*udAv,
	  // 	 dkdp(j, 0));
	}
      }

      for (size_t j = 0; j < nbasecells; j ++) {

	double udAvdT = 0.0;
	double udBvdT = 0.0;
	double udCvdT = 0.0;
	
	
	for (size_t i = 0; i < size; i ++) {
	  udAvdT += v(i, 0) * dAvdT(i, j);
	  udBvdT += v(i, 0) * dBvdT(i, j);
	  udCvdT += v(i, 0) * dCvdT(i, j);
	}

	dkdT(j, 0) = ((o2*udAvdT - udCvdT) - fundamental*udBvdT)/(2.0 * k * normB);
      }
	

      // fclose(fp);
    }

    return k;
  }

  real solve_fundamental_gradient_sep(const Mesh<real, maxorder> &mesh,
				      size_t boundaryorder,
				      real omega,
				      Spec1DMatrix<real> &dkdp,
				      Spec1DMatrix<real> &dUdp,
				      real &normA,
				      real &normB,
				      real &normC)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();
    Bs.setZero();
    
    real o2 = omega*omega;
    real scale = 0.0;
    
    for (size_t i = 0; i < size; i ++) {
      D(i, i) = o2*A(i, i);
      if (D(i, i) > scale) {
	scale = D(i, i);
      }
    }

    //scale = 1.0;
    
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) -= C(j, i);
	D(j, i) /= scale;
      }
      Bs(j, j) = B(j, j)/scale;
    }

    if (!SpecializedEigenProblem(D, Bs, work, eu, ev, lambda, Q, Z)) {
      ERROR("Failed to compute generalised eigen problem");
      return 0.0;
    }
    
    real fundamental = 0.0;
    v.resize(size, 1);
    u.resize(size, 1);

    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0 && k2 > fundamental) {
	  fundamental = k2;

	  for (size_t j = 0; j < size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    double unorm = 0.0;
    double vnorm = 0.0;
    for (size_t j = 0; j < size; j ++) {
      unorm += u(j, 0)*u(j, 0);
      vnorm += v(j, 0)*v(j, 0);
    }
    unorm = sqrt(unorm);
    vnorm = sqrt(vnorm);

    if (v(0, 0) > 0.0) {
      for (size_t j = 0; j < size; j ++) {
	u(j, 0) /= unorm;
	v(j, 0) /= vnorm;
      }
    } else {
      for (size_t j = 0; j < size; j ++) {
	u(j, 0) /= -unorm;
	v(j, 0) /= -vnorm;
      }
    }      
    
    real k = 0.0;
    
    if (fundamental > 0.0) {
      k = sqrt(fundamental);

      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      for (size_t j = 0; j < size; j ++) {
	normA += v(j, 0) * A(j, j) * v(j, 0);
	normB += v(j, 0) * B(j, j) * v(j, 0);
	real c = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  c += C(j, i) * v(i, 0);
	}
	normC += v(j, 0) * c;
      }

      size_t nbasecells;
      size_t nparameters = postcomputegradient(mesh,
					       boundaryorder,
					       v,
					       nbasecells);

      real galpha = normB/(2.0*k*omega*normA);
      real adv = 0.0;
      real bdv = 0.0;
      for (size_t i = 0; i < size; i ++) {
	adv += A(i, i) * v(i, 0);
	bdv += B(i, i) * v(i, 0);
      }

      //
      // Compute g_v (Gradient of U w.r.t. v)
      //
      gv.resize(size, 1);
      real t = 0.0;
      for (size_t i = 0; i < size; i ++) {
	gv(i, 0) = (2.0*k*(normA * B(i, i) - normB * A(i, i)) * v(i, 0))/(omega * normA * normA);
	t += gv(i, 0) * v(i, 0);
      }

      //
      // Compute component of g_v orthogonal to v
      //
      for (size_t i = 0; i < size; i ++) {
	gv(i, 0) -= t*v(i, 0);
      }

      //
      // Solve for lambda_0 with forward subsitution
      //
      adjointlambda0.resize(size, 1);
      SpecializedEigenProblemAdjoint<real>(D,
					   Bs,
					   Q,
					   Z,
					   fundamental, // k^2
					   gv,
					   adjointlambda0,
					   work);

      //
      // Ensure orthogonality
      //
      real lambda0du = 0.0;
      for (size_t i = 0; i < size; i ++) {
	lambda0du += adjointlambda0(i, 0) * u(i, 0);
      }
      for (size_t i = 0; i < size; i ++) {
	adjointlambda0(i, 0) -= lambda0du * u(i, 0);
      }
      

      adjointw.resize(size, 1);
      for (size_t i = 0; i < size; i ++) {
	adjointw(i, 0) = v(i, 0) * B(i, i)/normB;
      }

      real gamma = 0.0;
      for (size_t i = 0; i < size; i ++) {
	gamma -= adjointw(i, 0) * adjointlambda0(i, 0)/scale;
      }
      gamma -= galpha/normB;

      for (size_t i = 0; i < size; i ++) {
	adjointlambda0(i, 0) += gamma*u(i, 0)*scale;
      }

      //
      // Check (A^T - alpha B^T) lambda
      //
      double Bvlambda = 0.0;
      for (size_t i = 0; i < size; i ++) {
	double s = 0.0;
	for (size_t j = 0; j < size; j ++) {
	  s += (omega*omega*A(j, i) - C(j, i) - fundamental*B(i, j))/scale * adjointlambda0(j, 0);
	}

	Bvlambda += B(i, i) * v(i, 0) * adjointlambda0(i, 0);
      }
      
      //
      // -B v^T lambda
      //
      
      dkdp.resize(nparameters, 1);
      dkdp.setZero();

      dUdp.resize(nparameters, 1);
      dUdp.setZero();


      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;
	double ll = 0.0;
	
	for (size_t i = 0; i < size; i ++) {
	  udAv += v(i, 0) * dAv(i, j);
	  udBv += v(i, 0) * dBv(i, j);
	  udCv += v(i, 0) * dCv(i, j);

	  ll += adjointlambda0(i, 0) * (o2*dAv(i, j) - dCv(i, j) - fundamental*dBv(i, j))/scale;
	}

	dkdp(j, 0) = ((o2*udAv - udCv) - fundamental*udBv)/(2.0 * k * normB);

	dUdp(j, 0) = (normA*k*udBv - k*normB*udAv)/(omega*normA*normA) - ll;

      }
    }

    return k;
  }

  real solve_fundamental_gep_vector(real omega,
				    const Mesh<real, maxorder> &mesh,
				    size_t boundaryorder,
				    MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude,
				    real &L2norm)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i));
	Bs(j, i) = B(j, i);
      }
    }

    if (!GEP(D, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalised eigen problem");
      return 0.0;
    }
    
    real fundamental = 0.0;
    v.resize(size, 1);
    rA.resize(size, 1);
    rB.resize(size, 1);
    
    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0 && k2 > fundamental) {
	  fundamental = k2;

	  for (size_t j = 0; j < size; j ++) {
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    for (size_t i = 0; i < size; i ++) {
      real tA = 0.0;
      real tB = 0.0;

      for (size_t j = 0; j < size; j ++) {
	tA += (o2 * A(i, j) - C(i, j)) * v(j, 0);
	tB += B(i, j) * v(j, 0);
      }

      rA(i, 0) = tA;
      rB(i, 0) = tB;
    }

    fill_amplitude(mesh, v, boundaryorder, amplitude);

    L2norm = 0.0;
    real lastL2 = 0.0;
    real lastL = mesh.cells[0].nodes[0].L;
    for (size_t i = 0; i < amplitude.cells.size(); i ++) {

      real L2top = 0.0;
      real L2bottom = 0.0;
      for (size_t j = 0; j <= amplitude.cells[i].order; j ++) {
	L2top += Lobatto[amplitude.cells[i].order]->cardinal[j].dydx(-1.0)*amplitude.cells[i].nodes[j].V * 2.0/amplitude.cells[i].thickness;
	L2bottom += Lobatto[amplitude.cells[i].order]->cardinal[j].dydx(1.0)*amplitude.cells[i].nodes[j].V * 2.0/amplitude.cells[i].thickness;
      }

      L2top *= mesh.cells[i].nodes[0].L;
      L2bottom *= mesh.cells[i].nodes[mesh.cells[i].order].L;

      real meanL = (lastL + mesh.cells[i].nodes[mesh.cells[i].order].L)/2.0;
      real delta = (L2top - lastL2)/meanL;
      L2norm += (delta*delta);
      printf("%2d %16.9e %16.9e %16.9e\n", (int)i, L2top, L2bottom, delta);
      
      lastL2 = L2bottom;
      lastL = mesh.cells[i].nodes[mesh.cells[i].order].L;
    }

    if (mesh.boundary.rho > 0.0) {

      real L2 = 0.0;
      for (size_t j = 0; j <= amplitude.boundary.order; j ++) {
	L2 += Laguerre[amplitude.boundary.order]->cardinal[j].dydx(0.0)*amplitude.boundary.nodes[j].V * laguerrescale;
      }

      L2 *= mesh.boundary.L;
      real meanL = (lastL + mesh.boundary.L)/2.0;
      real delta = (L2 - lastL2)/meanL;
      L2norm += (delta*delta);
      printf("n  %16.9e %16.9e\n",
	     L2, delta);
    }

    L2norm = sqrt(L2norm);

    return sqrt(fundamental);
  }
  
  bool solve_all_gep(real omega, std::vector<real> &k)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i));
	Bs(j, i) = B(j, i);
      }
    }

    if (!GEP(D, Bs, work, eu, ev, lambda)) {
      FATAL("Failed to compute generalised eigen problem");
    }
    
    for (size_t i = 0; i < size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k2 = lambda(i, 0)/lambda(i, 2);

	if (k2 > 0.0) {

	  k.push_back(sqrt(k2));

	}
      }
    }

    return true;
  }

  void computeA(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    A.setZero();
    for (size_t i = 0; i < (ncells - 1); i ++) {

      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	A(offset + j, offset + j) +=
	  mesh.cells[i].thickness/2.0 *
	  mesh.cells[i].nodes[j].rho *
	  Lobatto[mesh.cells[i].order]->weights[j];

      }

      offset += mesh.cells[i].order;
    }

    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      // Do the last cell
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

      	A(offset + j, offset + j) +=
      	  mesh.cells[i].thickness/2.0 *
      	  mesh.cells[i].nodes[j].rho *
      	  Lobatto[mesh.cells[i].order]->weights[j];

      }
      offset += mesh.cells[i].order;

      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	A(offset + j, offset + j) +=
	  1.0/laguerrescale *
	  mesh.boundary.rho *
	  Laguerre[boundaryorder]->weights[j];
      }
      
      
      
    } else {
      //
      // Fixed
      //

      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

      	A(offset + j, offset + j) +=
      	  mesh.cells[i].thickness/2.0 *
      	  mesh.cells[i].nodes[j].rho *
      	  Lobatto[mesh.cells[i].order]->weights[j];

      }
      
    }
  }

  void computeAGradient(size_t nparameters,
			size_t nbasecells,
			const Mesh<real, maxorder> &mesh,
			size_t boundaryorder,
			Spec1DMatrix<real> &v)
  {
    dAv.resize(v.rows(), nparameters);
    dAv.setZero();

    dAvdT.resize(v.rows(), nbasecells);
    dAvdT.setZero();

    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    
    for (size_t i = 0; i < (ncells - 1); i ++) {
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	//
	// dA/dp = (thick)/2 * w * drho/dp
	//
	// jacobian = | drho0/dp0 dA0/dp0 ... dN0.dp0 drho1/dp0 ... |
	//            | drho0/dp1 dA0/dp1 ... dN0.dp1 drho1/dp1 ... |
	//            |    ...      ...                             |
	//

	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Col of dA/dp V, offset + j is row

	  dAv(offset + j, l) += mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] * 
	    mesh.cells[i].jacobian(k, j*6 + 0) *
	    v(offset + j, 0);

	}

	int l = mesh.cell_thickness_reference[i];
	dAvdT(offset + j, l) += 0.5 *
	   Lobatto[mesh.cells[i].order]->weights[j] *
	   v(offset + j, 0) *
	   mesh.cells[i].thickness_jacobian * 
	   mesh.cells[i].nodes[j].rho;
	
      }
      offset += mesh.cells[i].order;
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      // Do the last cell
      
      size_t i = mesh.cells.size() - 1;
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Row of dA/dp V, offset + j is col
	  
	  dAv(offset + j, l) += mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] * 
	    mesh.cells[i].jacobian(k, j*6 + 0) *
	    v(offset + j, 0);
	  
	}

	int l = mesh.cell_thickness_reference[i];
	dAvdT(offset + j, l) += 0.5 *
	   Lobatto[mesh.cells[i].order]->weights[j] *
	   v(offset + j, 0) *
	   mesh.cells[i].thickness_jacobian * 
	   mesh.cells[i].nodes[j].rho;
      }
      offset += mesh.cells[i].order;

      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {

	size_t nparametersincell = mesh.boundary_jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.boundary_parameter_offset + k; // Row of dA/dp V, offset + j is col

	  
	  dAv(offset + j, l) += 1.0/laguerrescale *
	    Laguerre[boundaryorder]->weights[j] *
	    mesh.boundary_jacobian(k, 0) * // Halfspace is constant so drho/dp stored only in first column
	    v(offset + j, 0);
	}
      }
      
    } else {
      //
      // Fixed
      //

      // Do the last cell
      size_t i = mesh.cells.size() - 1;
      size_t nparametersincell = mesh.cells[i].jacobian.rows();
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Row of dA/dp V, offset + j is col
	  
	  dAv(offset + j, l) += mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] * 
	    mesh.cells[i].jacobian(k, j*6 + 0) *
	    v(offset + j, 0);
	  
	}
      }
    }

    // printf("dAvdT\n");
    // for (int i = 0; i < v.rows(); i ++) {
    //   double t = 0.0;
    //   for (int j = 0; j < v.rows(); j ++) {
    // 	t += A(i, j) * v(j, 0);
    //   }
    //   printf("%16.9e %16.9e %16.9e\n", dAvdT(i, 0), dAvdT(i,0)/v(i,0), t/mesh.cells[0].thickness);
    // }
  }

  void computeB(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    B.setZero();

    for (size_t i = 0; i < (ncells - 1); i ++) {

      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	B(offset + j, offset + j) +=
	  mesh.cells[i].thickness/2.0 *
	  mesh.cells[i].nodes[j].N *
	  Lobatto[mesh.cells[i].order]->weights[j];

      }

      offset += mesh.cells[i].order;
    }

    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      // Do the last cell
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

      	B(offset + j, offset + j) +=
      	  mesh.cells[i].thickness/2.0 *
      	  mesh.cells[i].nodes[j].N *
      	  Lobatto[mesh.cells[i].order]->weights[j];

      }
      offset += mesh.cells[i].order;

      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	B(offset + j, offset + j) +=
	  1.0/laguerrescale *
	  mesh.boundary.N *
	  Laguerre[boundaryorder]->weights[j];
      }
      
    } else {
      //
      // Fixed
      //

      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	B(offset + j, offset + j) +=
	  mesh.cells[i].thickness/2.0 *
	  mesh.cells[i].nodes[j].N *
	  Lobatto[mesh.cells[i].order]->weights[j];

      }
      
    }
    
  }

  void computeBGradient(size_t nparameters,
			size_t nbasecells,
			const Mesh<real, maxorder> &mesh,
			size_t boundaryorder,
			Spec1DMatrix<real> &v)
  {
    dBv.resize(v.rows(), nparameters);
    dBv.setZero();

    dBvdT.resize(v.rows(), nbasecells);
    dBvdT.setZero();

    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    for (size_t i = 0; i < (ncells - 1); i ++) {
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	//
	// dB/dp = (thick)/2 * w * dN/dp
	//
	// jacobian = | drho0/dp0 dA0/dp0 ... dN0.dp0 drho1/dp0 ... |
	//            | drho0/dp1 dA0/dp1 ... dN0.dp1 drho1/dp1 ... |
	//            |    ...      ...                             |
	//
	
	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Row of dA/dp V, offset + j is col

	  dBv(offset + j, l) += mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] * 
	    mesh.cells[i].jacobian(k, j*6 + 5) *
	    v(offset + j, 0);

	}

	int l = mesh.cell_thickness_reference[i];
	dBvdT(offset + j, l) += 0.5 *
	   Lobatto[mesh.cells[i].order]->weights[j] *
	   v(offset + j, 0) *
	   mesh.cells[i].thickness_jacobian * 
	   mesh.cells[i].nodes[j].N;
	
      }
      offset += mesh.cells[i].order;
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      // Do the last cell
      size_t i = mesh.cells.size() - 1;
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Row of dA/dp V, offset + j is col
	  
	  dBv(offset + j, l) += mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] * 
	    mesh.cells[i].jacobian(k, j*6 + 5) *
	    v(offset + j, 0);
	  
	}

	int l = mesh.cell_thickness_reference[i];
	dBvdT(offset + j, l) += 0.5 *
	   Lobatto[mesh.cells[i].order]->weights[j] *
	   v(offset + j, 0) *
	   mesh.cells[i].thickness_jacobian * 
	   mesh.cells[i].nodes[j].N;
      }
      offset += mesh.cells[i].order;

      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {

	size_t nparametersincell = mesh.boundary_jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.boundary_parameter_offset + k; // Row of dA/dp V, offset + j is col

	  
	  dBv(offset + j, l) += 1.0/laguerrescale *
	    Laguerre[boundaryorder]->weights[j] *
	    mesh.boundary_jacobian(k, 5) * // Halfspace is constant so dN/dp stored only in fifth column
	    v(offset + j, 0);
	}
      }
      
    } else {
      //
      // Fixed
      //

      // Do the last cell
      size_t i = mesh.cells.size() - 1;
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Row of dA/dp V, offset + j is col
	  
	  dBv(offset + j, l) += mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] * 
	    mesh.cells[i].jacobian(k, j*6 + 5) *
	    v(offset + j, 0);
	  
	}
      }
    }

    // printf("dBvdT\n");
    // for (int i = 0; i < v.rows(); i ++) {
    //   double t = 0.0;
    //   for (int j = 0; j < v.rows(); j ++) {
    // 	t += B(i, j) * v(j, 0);
    //   }
    //   printf("%16.9e %16.9e %16.9e\n", dBvdT(i, 0), dBvdT(i,0)/v(i,0), t/mesh.cells[0].thickness);
    // }
  }

  void computeC(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    C.setZero();

    for (size_t i = 0; i < (ncells - 1); i ++) {

      for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	    C(offset + m, offset + n) += 
	      2.0/mesh.cells[i].thickness *
	      mesh.cells[i].nodes[j].L * 
	      Lobatto[mesh.cells[i].order]->weights[j] * 
	      Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
	      Lobatto[mesh.cells[i].order]->derivative_weights[n][j];

	  }
	}
      }
      offset += mesh.cells[i].order;
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      // Do the last cell
      size_t i = mesh.cells.size() - 1;
      
      for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	    
	    C(offset + m, offset + n) += 
	      2.0/mesh.cells[i].thickness *
	      mesh.cells[i].nodes[j].L * 
	      Lobatto[mesh.cells[i].order]->weights[j] * 
	      Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
	      Lobatto[mesh.cells[i].order]->derivative_weights[n][j];
	  }
	}
      }
      
      offset += mesh.cells[i].order;

      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t n = 0; n <= boundaryorder; n ++) {
	  for (size_t j = 0; j <= boundaryorder; j ++) {
	    
	    C(offset + m, offset + n) += 
	      laguerrescale *
	      mesh.boundary.L * 
	      Laguerre[boundaryorder]->weights[j] * 
	      Laguerre[boundaryorder]->derivative_weights[m][j] *
	      Laguerre[boundaryorder]->derivative_weights[n][j];
	    
	  }
	}
      }
      
    } else {
      //
      // Fixed
      //
      size_t i = ncells - 1;

      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	// Note subtle difference in use of < instead of <= as fixed boundary is zero displacement so
	// we remove it from the loop.
	for (size_t n = 0; n < mesh.cells[i].order; n ++) { 
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	    C(offset + m, offset + n) += 
	      2.0/mesh.cells[i].thickness *
	      mesh.cells[i].nodes[j].L * 
	      Lobatto[mesh.cells[i].order]->weights[j] * 
	      Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
	      Lobatto[mesh.cells[i].order]->derivative_weights[n][j];
	  }
	}
      }
    }
  }
  
  void computeCGradient(size_t nparameters,
			size_t nbasecells,
			const Mesh<real, maxorder> &mesh,
			size_t boundaryorder,
			Spec1DMatrix<real> &v)
  {
    dCv.resize(v.rows(), nparameters);
    dCv.setZero();

    dCvdT.resize(v.rows(), nbasecells);
    dCvdT.setZero();

    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    for (size_t i = 0; i < (ncells - 1); i ++) {

      for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;
	      dCv(offset + m, l) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].jacobian(k, j*6 + 4) * 
		Lobatto[mesh.cells[i].order]->weights[j] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		v(offset + n, 0);
	    }

	    int l = mesh.cell_thickness_reference[i];
	    dCvdT(offset + m, l) -= 
	      2.0/(mesh.cells[i].thickness * mesh.cells[i].thickness) *
	      mesh.cells[i].nodes[j].L * 
	      Lobatto[mesh.cells[i].order]->weights[j] * 
	      Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
	      Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
	      mesh.cells[i].thickness_jacobian *
	      v(offset + n, 0);
	    
	  }
	}
      }
      offset += mesh.cells[i].order;
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      // Do the last cell
      size_t i = mesh.cells.size() - 1;
      
      for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	    
	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;
	      dCv(offset + m, l) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].jacobian(k, j*6 + 4) * 
		Lobatto[mesh.cells[i].order]->weights[j] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		v(offset + n, 0);
	    }

	    int l = mesh.cell_thickness_reference[i];
	    dCvdT(offset + m, l) -= 
	      2.0/(mesh.cells[i].thickness * mesh.cells[i].thickness) *
	      mesh.cells[i].nodes[j].L * 
	      Lobatto[mesh.cells[i].order]->weights[j] * 
	      Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
	      Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
	      mesh.cells[i].thickness_jacobian *
	      v(offset + n, 0);
	    
	  }
	}
      }
      
      offset += mesh.cells[i].order;

      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t n = 0; n <= boundaryorder; n ++) {
	  for (size_t j = 0; j <= boundaryorder; j ++) {

	    size_t nparametersincell = mesh.boundary_jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.boundary_parameter_offset + k;
	      dCv(offset + m, l) += 
		laguerrescale *
		mesh.boundary_jacobian(k, 4) * 
		Laguerre[boundaryorder]->weights[j] * 
		Laguerre[boundaryorder]->derivative_weights[m][j] *
		Laguerre[boundaryorder]->derivative_weights[n][j] *
		v(offset + n, 0);
	    }
	  }
	}
      }
      
    } else {
      //
      // Fixed
      //
      size_t i = ncells - 1;

      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	// Note subtle difference in use of < instead of <= as fixed boundary is zero displacement so
	// we remove it from the loop.
	for (size_t n = 0; n < mesh.cells[i].order; n ++) { 
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;
	      dCv(offset + m, l) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].jacobian(k, j*6 + 4) * 
		Lobatto[mesh.cells[i].order]->weights[j] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		v(offset + n, 0);
	    }
	  }
	}
      }
    }

    // printf("dCvdT\n");
    // for (int i = 0; i < v.rows(); i ++) {
    //   double t = 0.0;
    //   for (int j = 0; j < v.rows(); j ++) {
    // 	t += C(i, j) * v(j, 0);
    //   }
    //   printf("%16.9e %16.9e %16.9e\n", dCvdT(i, 0), dCvdT(i,0)/v(i,0), -t/mesh.cells[0].thickness);
    // }
  }

  void computeD(real omega)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    D.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        D(j, i) = (o2 * A(j, i) - C(j, i))/B(j, j);
      }
    }
  }

  real computeDiagonalMatrixProductSum(Spec1DMatrix<real> &a,
				       Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < size; i ++) {
      s += a(i, i) * v(i, 0);
    }

    return s;
  }
  
  real computeSquareDiagonalMatrixProductSum(Spec1DMatrix<double> &u,
					     Spec1DMatrix<real> &a,
					     Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < size; i ++) {
      s += u(i, 0) * a(i, i) * v(i, 0);
    }

    return s;
  }

  real computeSquareMatrixProductSum(Spec1DMatrix<double> &u,
				     real omega,
				     Spec1DMatrix<real> &a,
				     Spec1DMatrix<real> &c,
				     Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < size; i ++) {
      real n = 0.0;
      for (size_t j = 0; j < size; j ++) {
	n += u(j, 0) * (omega*omega * a(j, i) - c(j, i));
      }
      
      s += n * v(i, 0);
    }

    return s;
  }

  real computeMatrixProductSum(Spec1DMatrix<real> &a,
			       Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < size; i ++) {
      for (size_t j = 0; j < size; j ++) {
	s += a(i, j) * v(i, 0);
      }
    }

    return s;
  }
  
  
  void printMatrix(FILE *fp, const Spec1DMatrix<real> &m) const
  {
    if (size > 0) {
      for (size_t j = 0; j < size; j ++) {
	for (size_t i = 0; i < size; i ++) {
	  
	  fprintf(fp, "%10.3g ", value(m(j, i)));
	}
	fprintf(fp, "\n");
      }
    }
  }
  
  void printA(FILE *fp) const
  {
    printMatrix(fp, A);
  }
  void printB(FILE *fp) const
  {
    printMatrix(fp, B);
  }
  void printC(FILE *fp) const
  {
    printMatrix(fp, C);
  }
  void printD(FILE *fp) const
  {
    printMatrix(fp, D);
  }

  void fill_amplitude(const Mesh<real, maxorder> &mesh,
		      const Spec1DMatrix<real> &eigen_vector,
		      size_t boundaryorder,
		      MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude)
  {
    int offset = 0;

    for (size_t i = 0; i < mesh.cells.size() - 1; i ++) {

      CellAmplitude<real, maxorder> ca;

      ca.thickness = mesh.cells[i].thickness;
      ca.order = mesh.cells[i].order;
      for (size_t j = 0; j <= ca.order; j ++) {
	ca.nodes[j].V = eigen_vector(offset + j, 0);
	ca.nodes[j].U = 0.0;
	ca.nodes[j].W = 0.0;
      }

      amplitude.cells.push_back(ca);

      offset += ca.order;
    }

    if (mesh.boundary.rho == 0.0) {

      //
      // Last cell
      //
      size_t i = mesh.cells.size() - 1;
      
      CellAmplitude<real, maxorder> ca;

      ca.thickness = mesh.cells[i].thickness;
      ca.order = mesh.cells[i].order;
      for (size_t j = 0; j < ca.order; j ++) {
	ca.nodes[j].V = eigen_vector(offset + j, 0);
	ca.nodes[j].U = 0.0;
	ca.nodes[j].W = 0.0;
      }
      ca.nodes[ca.order].V = 0.0; // Fixed boundary has fixed 0 at basement
      ca.nodes[ca.order].U = 0.0;
      ca.nodes[ca.order].W = 0.0;
      
      amplitude.cells.push_back(ca);
      
      //
      // Fixed boundary
      //
      amplitude.boundary.thickness = 0.0;
      amplitude.boundary.order = boundaryorder;
      for (size_t i = 0; i <= boundaryorder; i ++) {
	amplitude.boundary.nodes[i].V = 0.0;
	amplitude.boundary.nodes[i].U = 0.0;
	amplitude.boundary.nodes[i].W = 0.0;
      }

    } else {
      //
      // Last cell
      //
      size_t i = mesh.cells.size() - 1;
      
      CellAmplitude<real, maxorder> ca;

      ca.thickness = mesh.cells[i].thickness;
      ca.order = mesh.cells[i].order;
      
      for (size_t j = 0; j <= ca.order; j ++) {
	ca.nodes[j].V = eigen_vector(offset + j, 0);
	ca.nodes[j].U = 0.0;
	ca.nodes[j].W = 0.0;
      }

      amplitude.cells.push_back(ca);

      offset += ca.order;

      //
      // Laguerre halfspace
      //
      amplitude.boundary.thickness = 1.0;
      amplitude.boundary.order = boundaryorder;
      for (size_t j = 0; j <= boundaryorder; j ++) {
	amplitude.boundary.nodes[j].V = eigen_vector(offset + j, 0);
	amplitude.boundary.nodes[j].U = 0.0;
	amplitude.boundary.nodes[j].W = 0.0;
      }
    }
      
  }

  real interpolate_eigenvector(const MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude,
			       real z)
  {
    real dz = z;
    for (auto &c : amplitude.cells) {

      if (dz < c.thickness) {

	real xi = 2.0*dz/c.thickness - 1.0;

	if (xi < -1.0 || xi > 1.0) {
	  FATAL("xi out of range: %f", xi);
	}
	
	real V = 0.0;

	for (size_t i = 0; i <= c.order; i ++) {

	  V += Lobatto[c.order]->cardinal[i].value(xi)*c.nodes[i].V;

	}

	return V;
	
      } else {

	dz -= c.thickness;

      }
    }

    if (amplitude.boundary.thickness == 0.0) {
      //
      // Fixed boundary -> return 0
      //
      return 0.0;

    } else {
      //
      // Laguerre halfspace: dz = distance into halfspace, xi = dz*scale
      //
      real xi = dz * laguerrescale;

      real V = 0.0;
      for (size_t i = 0; i <= amplitude.boundary.order; i ++) {
	V += Laguerre[amplitude.boundary.order]->cardinal[i].value(xi)*amplitude.boundary.nodes[i].V;
      }

      return V;
    }
  }

  real interpolate_derived(const Mesh<real, maxorder> &mesh,
			   const MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude,
			   real z)
  {
    MeshParameter<real> p = mesh.interpolate(z);
    
    real dz = z;
    for (auto &c : amplitude.cells) {

      if (dz < c.thickness) {

	real xi = 2.0*dz/c.thickness - 1.0;

	if (xi < -1.0 || xi > 1.0) {
	  FATAL("xi out of range: %f", xi);
	}
	
	real V = 0.0;

	for (size_t i = 0; i <= c.order; i ++) {

	  V += Lobatto[c.order]->cardinal[i].dydx(xi)*c.nodes[i].V * 2.0/c.thickness;

	}

	return V * p.L;
	
      } else {

	dz -= c.thickness;

      }
    }

    if (amplitude.boundary.thickness == 0.0) {
      //
      // Fixed boundary -> return 0
      //
      return 0.0;

    } else {
      //
      // Laguerre halfspace: dz = distance into halfspace, xi = dz*scale
      //
      real xi = dz * laguerrescale;

      real V = 0.0;
      for (size_t i = 0; i <= amplitude.boundary.order; i ++) {
	V += Laguerre[amplitude.boundary.order]->cardinal[i].dydx(xi)*amplitude.boundary.nodes[i].V * laguerrescale;
      }

      return V * p.L;
    }
  }

  size_t size;

  std::array<LobattoQuadrature<double, maxorder>*, maxorder + 1> Lobatto;
  std::array<LaguerreQuadrature<double, maxboundaryorder>*, maxboundaryorder + 1> Laguerre;

  real laguerrescale;

  Spec1DMatrix<real> A;
  Spec1DMatrix<real> B;
  Spec1DMatrix<real> Bs;
  Spec1DMatrix<real> C;
  Spec1DMatrix<real> D;

  Spec1DMatrix<real> work;

  Spec1DMatrix<real> eu, ev;
  Spec1DMatrix<real> lambda;
  Spec1DMatrix<real> u;
  Spec1DMatrix<real> v;

  Spec1DMatrix<real> rA;
  Spec1DMatrix<real> rB;

  Spec1DMatrix<double> dB;
  Spec1DMatrix<double> dD;
  Spec1DMatrix<double> dwork;
  Spec1DMatrix<double> deu;
  Spec1DMatrix<double> dev;
  Spec1DMatrix<double> dlambda;
  Spec1DMatrix<double> du;
  Spec1DMatrix<double> dv;

  Spec1DMatrix<real> dAv;
  Spec1DMatrix<real> dBv;
  Spec1DMatrix<real> dCv;

  Spec1DMatrix<real> dAvdT;
  Spec1DMatrix<real> dBvdT;
  Spec1DMatrix<real> dCvdT;

  Spec1DMatrix<real> Q;
  Spec1DMatrix<real> Z;

  Spec1DMatrix<real> gv;
  Spec1DMatrix<real> adjointlambda0;
  Spec1DMatrix<real> adjointw;

  Spec1DMatrix<real> l2values;
};

#endif // lovematrices_hpp
  
