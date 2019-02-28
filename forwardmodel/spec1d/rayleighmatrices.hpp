#pragma once
#ifndef rayleighmatrices_hpp
#define rayleighmatrices_hpp

#include "mesh.hpp"

#include "lobattoquadrature.hpp"
#include "laguerrequadrature.hpp"

#include "spec1dmatrix.hpp"
#include "generalisedeigenproblem.hpp"
#include "generalsolve.hpp"
#include "specializedeigenproblem.hpp"

template
<
  typename real,
  size_t maxorder,
  size_t maxboundaryorder = maxorder
>
class RayleighMatrices {
public:

  RayleighMatrices() :
    size(0),
    IPIV(new int[1024]),
    IPIV_size(1024),
    disable_scale(false)
  {
    for (size_t i = 0; i <= maxorder; i ++) {
      Lobatto[i] = new LobattoQuadrature<double, maxorder>(i);
    }
    
    for (size_t i = 0; i <= maxboundaryorder; i ++) {
      Laguerre[i] = new LaguerreQuadrature<double, maxboundaryorder>(i);
    }
  }

  ~RayleighMatrices()
  {
    for (size_t i = 0; i <= maxorder; i ++) {
      delete Lobatto[i];
    }
    
    for (size_t i = 0; i <= maxboundaryorder; i ++) {
      delete Laguerre[i];
    }
  }

  void recompute(const Mesh<real, maxorder> &mesh, size_t boundaryorder, real scalex = 1.0, real scalez = 1.0)
  {
    if (boundaryorder > maxboundaryorder) {
      FATAL("Boundary order out of range: %d > %d", (int)boundaryorder, (int)maxboundaryorder);
    }

    size_t cell_nodes = mesh.nodes();

    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre boundary (extra maxorder - 1 nodes)
      //
      size = cell_nodes + boundaryorder;

    } else {

      //
      // Fixed boundary (last node is zero and removed from system of equations)
      //
      size = cell_nodes - 1;

    }

    laguerrescalex = scalex;
    laguerrescalez = scalez;

    A0x.resize(size, size);
    A0z.resize(size, size);
      
    Ax.resize(size, size);
    Bx.resize(size, size);
    Cx.resize(size, size);
    Dx.resize(size, size);

    Az.resize(size, size);
    Bz.resize(size, size);
    Cz.resize(size, size);
    Dz.resize(size, size);

    E.resize(4 * size, 4 * size);
    As.resize(4 * size, 4 * size);
    Bs.resize(4 * size, 4 * size);

    dAs.resize(4 * size, 4 * size);
    dBs.resize(4 * size, 4 * size);

    Ex.resize(size, size);
    Ez.resize(size, size);
    
    computeAx(mesh, boundaryorder);
    computeAz(mesh, boundaryorder);

    computeBx(mesh, boundaryorder);
    computeBz(mesh, boundaryorder);

    computeCx(mesh, boundaryorder);
    computeCz(mesh, boundaryorder);

    computeDx(mesh, boundaryorder);
    computeDz(mesh, boundaryorder);

  }

  size_t postcomputegradient(const Mesh<real, maxorder> &mesh,
			     size_t boundaryorder,
			     Spec1DMatrix<real> &v)
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
	
    computeAxGradient(nparameters, mesh, boundaryorder, v);
    computeAzGradient(nparameters, mesh, boundaryorder, v);
    computeBxGradient(nparameters, mesh, boundaryorder, v);
    computeBzGradient(nparameters, mesh, boundaryorder, v);
    computeCxGradient(nparameters, mesh, boundaryorder, v);
    computeCzGradient(nparameters, mesh, boundaryorder, v);
    computeDxGradient(nparameters, mesh, boundaryorder, v);
    computeDzGradient(nparameters, mesh, boundaryorder, v);

    return nparameters;
  }

  int get_size() const
  {
    return size;
  }
    

  real solve_fundamental_scaled(real omega,
				real &normA,
				real &normB,
				real &normC,
				real &normD,
				real &gamma,
				real &delta)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    gamma = computeE_scaled(omega, delta);
    
    if (!SpecializedEigenProblem(As, Bs, work, eu, ev, lambda, Q, Z)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }
    // if (!GEP(As, Bs, work, eu, ev, lambda)) {
    //   ERROR("Failed to compute generalized eigen problem");
    //   return 0.0;
    // }

    real sfundamental = 0.0;
    real fundamental = 0.0;
    v.resize(4*size, 1);
    u.resize(4*size, 1);

    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);
	if (absk > fundamental) {
	  fundamental = absk;
	  sfundamental = k;
	  for (size_t j = 0; j < 4*size; j ++) {
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
      normD = 0.0;
      
      for (size_t j = 0; j < size; j ++) {
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)

	normA +=
	  v(j, 0) * Ax(j, j) * v(j, 0) +
	  v(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  v(j, 0) * Bx(j, j) * v(j, 0) +
	  v(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += v(j, 0) * cx + v(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += v(j, 0) * dx + v(size + j, 0) * dz;
      }
    }
    
    return sfundamental * gamma;
  }

  real solve_fundamental_gradient(const Mesh<real, maxorder> &mesh,
				  size_t boundaryorder,
				  real omega,
				  Spec1DMatrix<real> &dkdp,
				  Spec1DMatrix<real> &dUdp,
				  real &normA,
				  real &normB,
				  real &normC,
				  real &normD)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);
    
    if (!GEP(As, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    real fundamental = 0.0;
    real sfundamental = 0.0;
    v.resize(4*size, 1);
    u.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);

	if (absk > fundamental) {
	  sfundamental = k;
	  fundamental = absk;

	  for (size_t j = 0; j < 4*size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    real k = 0.0;
    if (fundamental > 0.0) {
      size_t nparameters = postcomputegradient(mesh, boundaryorder, v);

      dkdp.resize(nparameters,1);
      dkdp.setZero();

      dUdp.resize(nparameters,1);
      dUdp.setZero();

      real norm = 0.0;
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      normD = 0.0;
      for (size_t j = 0; j < size; j ++) {
	norm -= u(j, 0) * (Bx(j, j) * gamma * gamma * delta) * v(j, 0) ;
	norm -= u(size + j, 0) * (Bz(j, j) * gamma * gamma * delta) * v(size + j, 0) ;
	norm -= u(2*size + j, 0) * v(2*size + j, 0);
	norm -= u(3*size + j, 0) * v(3*size + j, 0);

	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  u(j, 0) * Ax(j, j) * v(j, 0) +
	  u(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  u(j, 0) * Bx(j, j) * v(j, 0) +
	  u(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += u(j, 0) * cx + u(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += u(j, 0) * dx + u(size + j, 0) * dz;
      }

      k = sfundamental * gamma;

      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;
	double udDv = 0.0;
	
	for (size_t i = 0; i < size; i ++) {

	  dkdp(j, 0) += u(i, 0) * dCxv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(i, 0) * sfundamental * dBxv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += u(size + i, 0) * dCzv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(size + i, 0) * sfundamental * dBzv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += u(2*size + i, 0) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  dkdp(j, 0) += u(3*size + i, 0) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;

	  udAv += u(i, 0)*dAxv(i, j) + u(size + i, 0)*dAzv(i, j);
	  udBv += u(i, 0)*dBxv(i, j) + u(size + i, 0)*dBzv(i, j);
	  udCv += u(i, 0)*dCxv(i, j) + u(size + i, 0)*dCzv(i, j);
	  udDv += u(i, 0)*dDxv(i, j) + u(size + i, 0)*dDzv(i, j);
	  
	}

	dkdp(j, 0) *= gamma/norm;

	dUdp(j, 0) = (normA*(2.0*k*udBv + 2.0*dkdp(j, 0)*normB + udCv) - udAv*(2.0*k*normB + normC))/(2.0*omega*normA*normA);
      }

    }
    
    return k;
  }
  
  real solve_fundamental_gradient_v(const Mesh<real, maxorder> &mesh,
				    size_t boundaryorder,
				    real omega,
				    Spec1DMatrix<real> &dgxdv,
				    Spec1DMatrix<real> &dgzdv,
				    Spec1DMatrix<real> &dkdp,
				    Spec1DMatrix<real> &dUdp,
				    real &normA,
				    real &normB,
				    real &normC,
				    real &normD,
				    real &eH,
				    real &eV,
				    Spec1DMatrix<real> &dGvdp)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);

    adjointA.resize(4*size, 4*size);
    adjointB.resize(4*size, 2);
    adjointlambda.resize(4*size, 2);

    //
    // Precopy A matrix before it is destroyed by GEP
    //
    for (size_t i = 0;  i < 4*size; i ++) {
      for (size_t j = 0; j < 4*size; j ++) {
	adjointA(j, i) = As(i, j); // Note transpose
      }
    }
    
    if (!GEP(As, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    real fundamental = 0.0;
    real sfundamental = 0.0;
    v.resize(4*size, 1);
    u.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);

	if (absk > fundamental) {
	  sfundamental = k;
	  fundamental = absk;

	  for (size_t j = 0; j < 4*size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    real k = 0.0;
    if (fundamental > 0.0) {
      //
      // Normalize u, v
      //
      double unorm = 0.0;
      double vnorm = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	unorm += u(i, 0) * u(i, 0);
	vnorm += v(i, 0) * v(i, 0);
      }
      unorm = sqrt(unorm);
      vnorm = sqrt(vnorm);
      for (size_t i = 0; i < 4*size; i ++) {
	u(i, 0) /= unorm;
	v(i, 0) /= vnorm;
      }
	
      size_t nparameters = postcomputegradient(mesh, boundaryorder, v);

      dkdp.resize(nparameters,1);
      dkdp.setZero();

      dUdp.resize(nparameters,1);
      dUdp.setZero();

      real norm = 0.0;
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      normD = 0.0;
      for (size_t j = 0; j < size; j ++) {
	norm -= u(j, 0) * (Bx(j, j) * gamma * gamma * delta) * v(j, 0) ;
	norm -= u(size + j, 0) * (Bz(j, j) * gamma * gamma * delta) * v(size + j, 0) ;
	norm -= u(2*size + j, 0) * v(2*size + j, 0);
	norm -= u(3*size + j, 0) * v(3*size + j, 0);

	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  u(j, 0) * Ax(j, j) * v(j, 0) +
	  u(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  u(j, 0) * Bx(j, j) * v(j, 0) +
	  u(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += u(j, 0) * cx + u(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += u(j, 0) * dx + u(size + j, 0) * dz;
      }

      k = sfundamental * gamma;

      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;
	double udDv = 0.0;
	
	for (size_t i = 0; i < size; i ++) {

	  dkdp(j, 0) += u(i, 0) * dCxv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(i, 0) * sfundamental * dBxv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += u(size + i, 0) * dCzv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(size + i, 0) * sfundamental * dBzv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += u(2*size + i, 0) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  dkdp(j, 0) += u(3*size + i, 0) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;

	  udAv += u(i, 0)*dAxv(i, j) + u(size + i, 0)*dAzv(i, j);
	  udBv += u(i, 0)*dBxv(i, j) + u(size + i, 0)*dBzv(i, j);
	  udCv += u(i, 0)*dCxv(i, j) + u(size + i, 0)*dCzv(i, j);
	  udDv += u(i, 0)*dDxv(i, j) + u(size + i, 0)*dDzv(i, j);
	  
	}

	dkdp(j, 0) *= gamma/norm;

	dUdp(j, 0) = (normA*(2.0*k*udBv + 2.0*dkdp(j, 0)*normB + udCv) - udAv*(2.0*k*normB + normC))/(2.0*omega*normA*normA);
      }


      //
      // Adjoint matrix (initial set before call to GEP)
      //
      for (size_t j = 0; j < size; j ++) {
	adjointA(j, j) += Bx(j, j) * gamma * gamma * delta * sfundamental;
	adjointA(j + size, j + size) += Bz(j, j) * gamma * gamma * delta * sfundamental;
	adjointA(j + 2*size, j + 2*size) += sfundamental;
	adjointA(j + 3*size, j + 3*size) += sfundamental;
      }

      // printf("Adjoint A\n");
      // for (size_t j = 0; j < 4*size; j ++) {
      // 	printf("%02d %16.9e\n", adjointA(j, 0));
      // }
      
      //
      // Projection of gradient to orthogonal direction.
      //
      adjointB.setZero();
      // printf("gp\n");
      real vdgx = 0.0;
      real vdgz = 0.0;
      for (size_t i = 0; i < size; i ++) {
	vdgx += v(i, 0) * dgxdv(i, 0);
	vdgz += v(size + i, 0) * dgzdv(i, 0);
      }
      for (size_t i = 0; i < 4*size; i ++) {
	if (i < size) {
	  adjointB(i, 0) = dgxdv(i, 0);
	} else if (i < 2*size) {
	  adjointB(i, 1) = dgzdv(i - size, 0);
	}
	
	adjointB(i, 0) -= v(i, 0) * vdgx;
	adjointB(i, 1) -= v(i, 0) * vdgz;

	// printf("%02d %16.9e\n", adjointB(i, 0));
      }

      // v.T adjointB == 0
      // real ta = 0.0, tb = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	ta += adjointB(i, 0) * v(i, 0);
      // 	tb += adjointB(i, 1) * v(i, 0);
      // }
      
      //
      // Solve matrix equation
      //
      while (adjointA.rows() > IPIV_size) {
	delete IPIV;
	IPIV_size *= 2;
	IPIV = new int[IPIV_size];
      }
      adjointAt.resize(adjointA.rows(), adjointA.cols());
      
      for (size_t i = 0; i < 4*size; i ++) {
	for (size_t j = 0; j < 4*size; j ++) {
	  adjointAt(i, j) = adjointA(i, j);
	}
      }
      if (!GeneralSolve(adjointA, adjointB, IPIV)) {
	FATAL("Failed to solve adjoint");
      }

      adjointlambda.setZero();
      real uA = 0.0;
      real uB = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	uA += u(i, 0) * adjointB(i, 0);
	uB += u(i, 0) * adjointB(i, 1);
      }
      // printf("lambda\n");
      for (size_t i = 0; i < 4*size; i ++) {
	// adjointlambda(i, 0) = adjointB(i, 0);
	// adjointlambda(i, 1) = adjointB(i, 1);
	adjointlambda(i, 0) = adjointB(i, 0) - u(i, 0)*uA;
	adjointlambda(i, 1) = adjointB(i, 1) - u(i, 0)*uB;
	// printf("%02d %16.9e\n", i, adjointlambda(i, 0));
      }

      //
      // Adjoint lambda is now lambda_0, need to add correction for solution
      //
      //
      // w = vT B (B diagonal)
      //
      real wnorm = 0.0;
      adjointw.resize(4*size, 1);
      for (size_t i = 0; i < size; i ++) {

	adjointw(i, 0) = v(i, 0) * Bx(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i, 0)*adjointw(i, 0);

	adjointw(i + size, 0) = v(i + size, 0) * Bz(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i + size, 0)*adjointw(i + size, 0);

	adjointw(i + 2*size, 0) = v(i + 2*size, 0);
	wnorm += adjointw(i + 2*size, 0)*adjointw(i + 2*size, 0);

	adjointw(i + 3*size, 0) = v(i + 3*size, 0);
	wnorm += adjointw(i + 3*size, 0)*adjointw(i + 3*size, 0);

      }

      real dA_a = 0.0;
      real dA_b = 0.0;
      real dB = 0.0;
      wnorm = sqrt(wnorm);
      
      for (size_t i = 0; i < 4*size; i ++) {
	adjointw(i, 0) /= wnorm;
	dA_a += adjointw(i, 0) * adjointlambda(i, 0);
	dA_b += adjointw(i, 0) * adjointlambda(i, 1);
	dB += adjointw(i, 0) * u(i, 0);
      }

      if (dB == 0.0) {
	FATAL("Adjoint w and u are orthogonal");
      }
      real adjointscale_a = -dA_a/dB;
      real adjointscale_b = -dA_b/dB;

      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) += adjointscale_a * u(i, 0);
	adjointlambda(i, 1) += adjointscale_b * u(i, 0);
      }

      // printf("A * lambda\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real c = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  c += adjointAt(i, j) * adjointlambda(j, 0);
      // 	}
      // 	printf("%02d %16.9e\n", i, c);
      // }
      
      // v.T adjointlambda == 0
      // real sa = 0.0, sb = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	sa += adjointlambda(i, 0) * v(i, 0);
      // 	sb += adjointlambda(i, 1) * v(i, 0);
      // }
      // printf("t = %e, %e, s = %e, %e\n", ta, tb, sa, sb);
      
      dGvdp.resize(nparameters, 2);

      for (size_t k = 0; k < 2; k ++) {
	for (size_t j = 0; j < nparameters; j ++) {

	  
	  double a = 0.0;
	  for (size_t i = 0; i < size; i ++) {
	    
	    a += adjointlambda(i, k)*dCxv(i, j) * delta * gamma;
	    a += adjointlambda(size + i, k) * dCzv(i, j) * delta * gamma;
	    // a += adjointlambda(2*size + i, 0) * (omega*omega*dAxv(i, j) - dDxv(i, j)) * delta;
	    // a += adjointlambda(3*size + i, 0) * (omega*omega*dAzv(i, j) - dDzv(i, j)) * delta;

	    a += adjointlambda(2*size + i, k) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	    a += adjointlambda(3*size + i, k) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;
	    
	    a += sfundamental*adjointlambda(i, k)*dBxv(i, j) * delta * gamma * gamma;
	    a += sfundamental*adjointlambda(size + i, k)*dBzv(i, j) * delta * gamma * gamma;
	  }

	  dGvdp(j, k) = -a;
	}
      }
      
      eH = v(0, 0);
      eV = v(size, 0);
    }
    
    return k;
  }

  real solve_fundamental_gradient_group(const Mesh<real, maxorder> &mesh,
					size_t boundaryorder,
					real omega,
					Spec1DMatrix<real> &dkdp,
					Spec1DMatrix<real> &dUdp,
					real &normA,
					real &normB,
					real &normC,
					real &normD)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);

    adjointA.resize(4*size, 4*size);
    adjointB.resize(4*size, 1);
    adjointlambda.resize(4*size, 1);

    //
    // Precopy A matrix before it is destroyed by GEP
    //
    // for (size_t i = 0;  i < 4*size; i ++) {
    //   for (size_t j = 0; j < 4*size; j ++) {
    // 	adjointA(j, i) = As(i, j); // Note transpose
    //   }
    // }
    
    if (!SpecializedEigenProblem(As, Bs, work, eu, ev, lambda, Q, Z)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    real fundamental = 0.0;
    real sfundamental = 0.0;
    v.resize(4*size, 1);
    u.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);

	if (k > 0.0 && absk > fundamental) {
	  sfundamental = k;
	  fundamental = absk;

	  for (size_t j = 0; j < 4*size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    real k = 0.0;
    if (fundamental > 0.0) {
      //
      // Normalize u, v
      //
      double unorm = 0.0;
      double vnorm = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	unorm += u(i, 0) * u(i, 0);
	vnorm += v(i, 0) * v(i, 0);
      }
      unorm = sqrt(unorm);
      vnorm = sqrt(vnorm);
      for (size_t i = 0; i < 4*size; i ++) {
	u(i, 0) /= unorm;
	v(i, 0) /= vnorm;
      }
	
      size_t nparameters = postcomputegradient(mesh, boundaryorder, v);

      k = sfundamental * gamma;

      dkdp.resize(nparameters,1);
      dkdp.setZero();

      dUdp.resize(nparameters,1);
      dUdp.setZero();

      real norm = 0.0;
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      normD = 0.0;
      for (size_t j = 0; j < size; j ++) {
	norm -= u(j, 0) * (Bx(j, j) * gamma * gamma * delta) * v(j, 0) ;
	norm -= u(size + j, 0) * (Bz(j, j) * gamma * gamma * delta) * v(size + j, 0) ;
	norm -= u(2*size + j, 0) * v(2*size + j, 0);
	norm -= u(3*size + j, 0) * v(3*size + j, 0);

	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  v(j, 0) * Ax(j, j) * v(j, 0) +
	  v(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  v(j, 0) * Bx(j, j) * v(j, 0) +
	  v(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += v(j, 0) * cx + v(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += v(j, 0) * dx + v(size + j, 0) * dz;
      }

      real galpha = (gamma * normB)/(omega * normA);
      


      //
      // Adjoint matrix (initial set before call to GEP)
      //
      // for (size_t j = 0; j < size; j ++) {
      // 	adjointA(j, j) += Bx(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + size, j + size) += Bz(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + 2*size, j + 2*size) += sfundamental;
      // 	adjointA(j + 3*size, j + 3*size) += sfundamental;
      // }

      // sepdump("debugAs.txt", adjointA);
      
      // printf("Adjoint A\n");
      // for (size_t j = 0; j < 4*size; j ++) {
      // 	printf("%02d %16.9e\n", adjointA(j, 0));
      // }

      //
      // Set g_v to (1/2omega) (Anorm*CV - CnormAV + 2k Bnorm AV - 2k Anorm BV)/Anorm^2
      //
      adjointB.setZero();
      real t = 0.0;
      for (size_t i = 0; i < size; i ++) {
	real av1 = 0.0;
	real av2 = 0.0;
	real bv1 = 0.0;
	real bv2 = 0.0;
	real cv1 = 0.0;
	real cv2 = 0.0;
	  
	for (size_t j = 0; j < size; j ++) {
	  av1 += Ax(i, j) * v(j, 0);
	  av2 += Az(i, j) * v(size + j, 0);
	  bv1 += Bx(i, j) * v(j, 0);
	  bv2 += Bz(i, j) * v(size + j, 0);
	  cv1 += Cx(i, j) * v(size + j, 0);
	  cv2 += Cz(i, j) * v(j, 0);
	}
	
	adjointB(i, 0) +=
	  (normA*cv1 - normC*av1 - 2.0*k*normB*av1 + 2.0*k*normA*bv1)/
	  (normA * normA * 2.0 * omega);
	
	adjointB(size + i, 0) +=
	  (normA*cv2 - normC*av2 - 2.0*k*normB*av2 + 2.0*k*normA*bv2)/
	  (normA * normA * 2.0 * omega);

	t +=
	  adjointB(i, 0) * v(i, 0) +
	  adjointB(size + i, 0) * v(size + i, 0);
      }

      // printf("adjointB\n");
      // for (size_t i = 0; i < 2*size; i ++) {
      // 	printf("%16.9e\n", adjointB(i, 0));
      // }
	       
      //
      // Projection of gradient to orthogonal direction.
      //
      for (size_t i = 0; i < 4*size; i ++) {
      	adjointB(i, 0) -= t * v(i, 0);
      }
      
      //
      // Solve for lambda_0 with forward substitution
      //
      adjointlambda.resize(4*size, 1);
      if (!SpecializedEigenProblemAdjoint<real>(As,
						Bs,
						Q,
						Z,
						sfundamental,
						adjointB,
						adjointlambda,
						work)) {
	ERROR("Failed to forward subsitute");
	return 0.0;
      }

      // printf("Lambda0\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	printf("%16.9e %16.9e\n", adjointlambda(i, 0), adjointlambda(i, 1));
      // }
      
      //
      // Check
      //
      // printf("Check a1 == b2 : a2 == b2 (using forward subs)\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real a = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  a += adjointA(i, j) * adjointlambda(j, 0);
      // 	}
      // 	printf("%02d %16.9e %16.9e\n",
      // 	       (int)i,
      // 	       a, adjointB(i, 0));
      // }

      //
      // Check
      //
      // printf("Check 2\n");
      // while (adjointA.rows() > IPIV_size) {
      // 	delete IPIV;
      // 	IPIV_size *= 2;
      // 	IPIV = new int[IPIV_size];
      // }
      // adjointAt.resize(adjointA.rows(), adjointA.cols());
      // adjointBt.resize(adjointB.rows(), adjointB.cols());
      
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  adjointAt(i, j) = adjointA(i, j);
      // 	}
      // }
      
      // 	for (size_t j = 0; j < adjointB.cols(); j ++) {
      // 	  adjointBt(i, j) = adjointB(i, j);
      // 	}
      // }
      
      // if (!GeneralSolve(adjointA, adjointB, IPIV)) {
      // 	FATAL("Failed to solve adjoint");
      // }
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	printf("%2d %16.9e %16.9e : %16.9e %16.9e\n",
      // 	       (int)i,
      // 	       adjointlambda(i, 0), adjointB(i, 0),
      // 	       adjointlambda(i, 1), adjointB(i, 1));
      // }

      //
      // Check 3
      //
      // printf("Check 3: a1 == b1 : a2 == b2 (using inverse)\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real a = 0.0;
      // 	real b = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  a += adjointAt(i, j) * adjointB(j, 0);
      // 	  b += adjointAt(i, j) * adjointB(j, 1);
      // 	}
      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       a, adjointBt(i, 0),
      // 	       b, adjointBt(i, 1));
      // }
      
      //
      // Ensure orthogonality
      //
      real uA = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	uA += u(i, 0) * adjointlambda(i, 0);
      }
      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) -= u(i, 0)*uA;
      }

      // printf("Check 4: a1 == b1 : a2 == b2 projected\n");
      // uA = 0.0;
      // uB = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	uA += u(i, 0) * adjointB(i, 0);
      // 	uB += u(i, 0) * adjointB(i, 1);
      // }
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	adjointB(i, 0) -= u(i, 0)*uA;
      // 	adjointB(i, 1) -= u(i, 0)*uB;

      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       adjointlambda(i, 0), adjointB(i, 0),
      // 	       adjointlambda(i, 1), adjointB(i, 1));
      // }
      
      

      //
      // Adjoint lambda is now lambda_0, need to add correction for solution
      //
      //
      // w = vT B (B diagonal)
      //
      real wnorm = 0.0;
      adjointw.resize(4*size, 1);
      for (size_t i = 0; i < size; i ++) {

	adjointw(i, 0) = v(i, 0) * Bx(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i, 0)*adjointw(i, 0);
	
	adjointw(i + size, 0) = v(i + size, 0) * Bz(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i + size, 0)*adjointw(i + size, 0);

	adjointw(i + 2*size, 0) = v(i + 2*size, 0);
	wnorm += adjointw(i + 2*size, 0)*adjointw(i + 2*size, 0);
	adjointw(i + 3*size, 0) = v(i + 3*size, 0);
	wnorm += adjointw(i + 3*size, 0)*adjointw(i + 3*size, 0);
      }

      wnorm = sqrt(wnorm);
      real beta = 0.0;
      real wdotu = 0.0;
      
      for (size_t i = 0; i < 4*size; i ++) {
	adjointw(i, 0) /= wnorm;
	beta += adjointw(i, 0) * adjointlambda(i, 0);
	wdotu += adjointw(i, 0) * u(i, 0);
      }

      if (wdotu == 0.0) {
	FATAL("Adjoint w and u are orthogonal");
      }

      beta = (-beta - galpha/wnorm)/wdotu;
      
      // printf("beta %16.9e (%16.9e %16.9e)\n", beta, galpha, normB);
      
      // printf("Check A*lambda0 = (1 - vvT)g_v\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      //  	real c = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      //  	  c += adjointAt(i, j) * adjointlambda(j, 0);
      //  	}
      //  	printf("%02d %16.9e %16.9e \n", (int)i, c, adjointB(i, 0));
      // }

      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) += beta * u(i, 0);
      }

      // printf("Check A*lambda = (1 - vvT)g_v\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      //  	real c = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      //  	  c += adjointAt(i, j) * adjointlambda(j, 0);
      //  	}
      //  	printf("%02d %16.9e %16.9e (%16.9e %16.9e %16.9e %16.9e)\n", (int)i, c, adjointB(i, 0), v(i, 0), u(i, 0), adjointw(i, 0), adjointlambda(i, 0));
      // }
      
      // v.T adjointlambda == 0
      // real sa = 0.0, sb = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	sa += adjointlambda(i, 0) * v(i, 0);
      // 	sb += adjointlambda(i, 1) * v(i, 0);
      // }
      // printf("t = %e, %e, s = %e, %e\n", ta, tb, sa, sb);

      //
      // Check -(B v)^T lambda = g_alpha
      //
      real bvTlambda = 0.0;
      for (size_t i = 0; i < size; i ++) {
      	bvTlambda += gamma*gamma*delta*Bx(i, i) * v(i, 0) * adjointlambda(i, 0);
      	bvTlambda += gamma*gamma*delta*Bz(i, i) * v(size + i, 0) * adjointlambda(size + i, 0);
      	bvTlambda += v(2*size + i, 0) * adjointlambda(2*size + i, 0);
      	bvTlambda += v(3*size + i, 0) * adjointlambda(3*size + i, 0);
      }
      printf("Check BvT Lambda: %16.9e %16.9e %16.9e\n", bvTlambda, galpha, sfundamental);
      
      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;
	double udDv = 0.0;
	double a = 0.0;
	
	for (size_t i = 0; i < size; i ++) {
	  
	  dkdp(j, 0) += u(i, 0) * dCxv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(i, 0) * sfundamental * dBxv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += u(size + i, 0) * dCzv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(size + i, 0) * sfundamental * dBzv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += u(2*size + i, 0) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  dkdp(j, 0) += u(3*size + i, 0) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;

	  udAv += v(i, 0)*dAxv(i, j) + v(size + i, 0)*dAzv(i, j);
	  udBv += v(i, 0)*dBxv(i, j) + v(size + i, 0)*dBzv(i, j);
	  udCv += v(i, 0)*dCxv(i, j) + v(size + i, 0)*dCzv(i, j);
	  udDv += v(i, 0)*dDxv(i, j) + v(size + i, 0)*dDzv(i, j);

	  a += adjointlambda(i, 0)*dCxv(i, j) * delta * gamma;
	  a += adjointlambda(size + i, 0) * dCzv(i, j) * delta * gamma;
	  // a += adjointlambda(2*size + i, 0) * (omega*omega*dAxv(i, j) - dDxv(i, j)) * delta;
	  // a += adjointlambda(3*size + i, 0) * (omega*omega*dAzv(i, j) - dDzv(i, j)) * delta;
	  
	  a += adjointlambda(2*size + i, 0) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  a += adjointlambda(3*size + i, 0) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;
	  
	  a += sfundamental*adjointlambda(i, 0)*dBxv(i, j) * delta * gamma * gamma;
	  a += sfundamental*adjointlambda(size + i, 0)*dBzv(i, j) * delta * gamma * gamma;
	}

	dkdp(j, 0) *= gamma/norm;

	// Need to add in partials for matrix derivatives here
	// Eg from love : (normA*k*udBv - k*normB*udAv)/(omega*normA*normA)
	
	dUdp(j, 0) =
	  (normA*(udCv + 2.0*k*udBv) - udAv*(normC + 2.0*k*normB))/(2.0*omega*normA*normA) -
	  a;
      }
    }
    
    return k;
  }

  real solve_fundamental_gradient_v2(const Mesh<real, maxorder> &mesh,
				    size_t boundaryorder,
				    real omega,
				    Spec1DMatrix<real> &dgxdv,
				    Spec1DMatrix<real> &dgzdv,
				    Spec1DMatrix<real> &dkdp,
				    Spec1DMatrix<real> &dUdp,
				    real &normA,
				    real &normB,
				    real &normC,
				    real &normD,
				    real &eH,
				    real &eV,
				    Spec1DMatrix<real> &dGvdp)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);

    adjointA.resize(4*size, 4*size);
    adjointB.resize(4*size, 2);
    adjointlambda.resize(4*size, 2);

    //
    // Precopy A matrix before it is destroyed by GEP
    //
    // for (size_t i = 0;  i < 4*size; i ++) {
    //   for (size_t j = 0; j < 4*size; j ++) {
    // 	adjointA(j, i) = As(i, j); // Note transpose
    //   }
    // }
    
    if (!SpecializedEigenProblem(As, Bs, work, eu, ev, lambda, Q, Z)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    real fundamental = 0.0;
    real sfundamental = 0.0;
    v.resize(4*size, 1);
    u.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);

	if (absk > fundamental) {
	  sfundamental = k;
	  fundamental = absk;

	  for (size_t j = 0; j < 4*size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    real k = 0.0;
    if (fundamental > 0.0) {
      //
      // Normalize u, v
      //
      double unorm = 0.0;
      double vnorm = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	unorm += u(i, 0) * u(i, 0);
	vnorm += v(i, 0) * v(i, 0);
      }
      unorm = sqrt(unorm);
      vnorm = sqrt(vnorm);
      for (size_t i = 0; i < 4*size; i ++) {
	u(i, 0) /= unorm;
	v(i, 0) /= vnorm;
      }
	
      size_t nparameters = postcomputegradient(mesh, boundaryorder, v);

      dkdp.resize(nparameters,1);
      dkdp.setZero();

      dUdp.resize(nparameters,1);
      dUdp.setZero();

      real norm = 0.0;
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      normD = 0.0;
      for (size_t j = 0; j < size; j ++) {
	norm -= u(j, 0) * (Bx(j, j) * gamma * gamma * delta) * v(j, 0) ;
	norm -= u(size + j, 0) * (Bz(j, j) * gamma * gamma * delta) * v(size + j, 0) ;
	norm -= u(2*size + j, 0) * v(2*size + j, 0);
	norm -= u(3*size + j, 0) * v(3*size + j, 0);

	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  v(j, 0) * Ax(j, j) * v(j, 0) +
	  v(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  v(j, 0) * Bx(j, j) * v(j, 0) +
	  v(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += v(j, 0) * cx + v(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += v(j, 0) * dx + v(size + j, 0) * dz;
      }

      k = sfundamental * gamma;

      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;
	double udDv = 0.0;
	
	for (size_t i = 0; i < size; i ++) {

	  dkdp(j, 0) += v(i, 0) * dCxv(i, j) * delta * gamma;
	  dkdp(j, 0) += v(i, 0) * sfundamental * dBxv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += v(size + i, 0) * dCzv(i, j) * delta * gamma;
	  dkdp(j, 0) += v(size + i, 0) * sfundamental * dBzv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += v(2*size + i, 0) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  dkdp(j, 0) += v(3*size + i, 0) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;

	  udAv += u(i, 0)*dAxv(i, j) + u(size + i, 0)*dAzv(i, j);
	  udBv += u(i, 0)*dBxv(i, j) + u(size + i, 0)*dBzv(i, j);
	  udCv += u(i, 0)*dCxv(i, j) + u(size + i, 0)*dCzv(i, j);
	  udDv += u(i, 0)*dDxv(i, j) + u(size + i, 0)*dDzv(i, j);
	  
	}

	dkdp(j, 0) *= gamma/norm;

	dUdp(j, 0) = (normA*(2.0*k*udBv + 2.0*dkdp(j, 0)*normB + udCv) - udAv*(2.0*k*normB + normC))/(2.0*omega*normA*normA);
      }


      //
      // Adjoint matrix (initial set before call to GEP)
      //
      // for (size_t j = 0; j < size; j ++) {
      // 	adjointA(j, j) += Bx(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + size, j + size) += Bz(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + 2*size, j + 2*size) += sfundamental;
      // 	adjointA(j + 3*size, j + 3*size) += sfundamental;
      // }

      // sepdump("debugAs.txt", adjointA);
      
      // printf("Adjoint A\n");
      // for (size_t j = 0; j < 4*size; j ++) {
      // 	printf("%02d %16.9e\n", adjointA(j, 0));
      // }
      
      //
      // Projection of gradient to orthogonal direction.
      //
      adjointB.setZero();
      // printf("gp\n");
      real vdgx = 0.0;
      real vdgz = 0.0;
      for (size_t i = 0; i < size; i ++) {
	vdgx += v(i, 0) * dgxdv(i, 0);
	vdgz += v(size + i, 0) * dgzdv(i, 0);
      }
      for (size_t i = 0; i < 4*size; i ++) {
	if (i < size) {
	  adjointB(i, 0) = dgxdv(i, 0);
	} else if (i < 2*size) {
	  adjointB(i, 1) = dgzdv(i - size, 0);
	}
	
	adjointB(i, 0) -= v(i, 0) * vdgx;
	adjointB(i, 1) -= v(i, 0) * vdgz;

	// printf("%02d %16.9e\n", adjointB(i, 0));
      }

      // v.T adjointB == 0
      // real ta = 0.0, tb = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	ta += adjointB(i, 0) * v(i, 0);
      // 	tb += adjointB(i, 1) * v(i, 0);
      // }
      // SpecializedEigenProblemVerify(As, Q, Z, tAs, work);
      // SpecializedEigenProblemVerify(Bs, Q, Z, tBs, work);

      // printf("Check back\n");
      // double maxdiff = 0.0;
      // for (int i = 0; i < 4*size; i ++) {
      // 	for (int j = 0; j < 4*size; j ++) {

      // 	  double a = tAs(i, j) - sfundamental*tBs(i, j);
      // 	  double b = adjointA(j, i);

      // 	  double diff = fabs(a - b);
      // 	  if (diff > maxdiff) {
      // 	    maxdiff = diff;
      // 	  }
      // 	  if (diff > 1.0) {
      // 	    printf("  %2d, %2d : %16.9e\n", i, j, diff);
      // 	  }
      // 	}
      // }
      // printf("  max %16.9e\n", maxdiff );

      //
      // Solve for lambda_0 with forward substitution
      //
      adjointlambda.resize(4*size, 2);
      if (!SpecializedEigenProblemAdjoint<real>(As,
						Bs,
						Q,
						Z,
						sfundamental,
						adjointB,
						adjointlambda,
						work)) {
	ERROR("Failed to forward subsitute");
	return 0.0;
      }

      // printf("Lambda0\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	printf("%16.9e %16.9e\n", adjointlambda(i, 0), adjointlambda(i, 1));
      // }
      
      //
      // Check
      //
      // printf("Check a1 == b2 : a2 == b2 (using forward subs)\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real a = 0.0;
      // 	real b = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  a += adjointA(i, j) * adjointlambda(j, 0);
      // 	  b += adjointA(i, j) * adjointlambda(j, 1);
      // 	}
      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       a, adjointB(i, 0),
      // 	       b, adjointB(i, 1));
      // }

      //
      // Check
      //
      // printf("Check 2\n");
      // while (adjointA.rows() > IPIV_size) {
      // 	delete IPIV;
      // 	IPIV_size *= 2;
      // 	IPIV = new int[IPIV_size];
      // }
      // adjointAt.resize(adjointA.rows(), adjointA.cols());
      // adjointBt.resize(adjointB.rows(), adjointB.cols());
      
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  adjointAt(i, j) = adjointA(i, j);
      // 	}

      // 	for (size_t j = 0; j < adjointB.cols(); j ++) {
      // 	  adjointBt(i, j) = adjointB(i, j);
      // 	}
      // }
      
      // if (!GeneralSolve(adjointA, adjointB, IPIV)) {
      // 	FATAL("Failed to solve adjoint");
      // }
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	printf("%2d %16.9e %16.9e : %16.9e %16.9e\n",
      // 	       (int)i,
      // 	       adjointlambda(i, 0), adjointB(i, 0),
      // 	       adjointlambda(i, 1), adjointB(i, 1));
      // }

      //
      // Check 3
      //
      // printf("Check 3: a1 == b1 : a2 == b2 (using inverse)\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real a = 0.0;
      // 	real b = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  a += adjointAt(i, j) * adjointB(j, 0);
      // 	  b += adjointAt(i, j) * adjointB(j, 1);
      // 	}
      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       a, adjointBt(i, 0),
      // 	       b, adjointBt(i, 1));
      // }
      
      //
      // Ensure lambda0 orthogonal (Need to be orthogonal to u as non-symmetric)
      //
      real uA = 0.0;
      real uB = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	uA += u(i, 0) * adjointlambda(i, 0);
	uB += u(i, 0) * adjointlambda(i, 1);
      }
      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) -= u(i, 0)*uA;
	adjointlambda(i, 1) -= u(i, 0)*uB;
      }

      // printf("Check 4: a1 == b1 : a2 == b2 projected\n");
      // uA = 0.0;
      // uB = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	uA += u(i, 0) * adjointB(i, 0);
      // 	uB += u(i, 0) * adjointB(i, 1);
      // }
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	adjointB(i, 0) -= u(i, 0)*uA;
      // 	adjointB(i, 1) -= u(i, 0)*uB;

      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       adjointlambda(i, 0), adjointB(i, 0),
      // 	       adjointlambda(i, 1), adjointB(i, 1));
      // }
      
      

      //
      // Adjoint lambda is now lambda_0, need to add correction for solution
      //
      //
      // w = vT B (B diagonal)
      //
      real wnorm = 0.0;
      adjointw.resize(4*size, 1);
      for (size_t i = 0; i < size; i ++) {

	adjointw(i, 0) = v(i, 0) * Bx(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i, 0)*adjointw(i, 0);

	adjointw(i + size, 0) = v(i + size, 0) * Bz(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i + size, 0)*adjointw(i + size, 0);

	adjointw(i + 2*size, 0) = v(i + 2*size, 0);
	wnorm += adjointw(i + 2*size, 0)*adjointw(i + 2*size, 0);

	adjointw(i + 3*size, 0) = v(i + 3*size, 0);
	wnorm += adjointw(i + 3*size, 0)*adjointw(i + 3*size, 0);

      }

      real dA_a = 0.0;
      real dA_b = 0.0;
      real dB = 0.0;
      wnorm = sqrt(wnorm);
      
      for (size_t i = 0; i < 4*size; i ++) {
	adjointw(i, 0) /= wnorm;
	dA_a += adjointw(i, 0) * adjointlambda(i, 0);
	dA_b += adjointw(i, 0) * adjointlambda(i, 1);
	dB += adjointw(i, 0) * u(i, 0);
      }

      if (dB == 0.0) {
	FATAL("Adjoint w and u are orthogonal");
      }
      real adjointscale_a = -dA_a/dB;
      real adjointscale_b = -dA_b/dB;

      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) += adjointscale_a * u(i, 0);
	adjointlambda(i, 1) += adjointscale_b * u(i, 0);
      }

      // printf("A * lambda\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real c = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  c += adjointAt(i, j) * adjointlambda(j, 0);
      // 	}
      // 	printf("%02d %16.9e\n", i, c);
      // }
      
      // v.T adjointlambda == 0
      // real sa = 0.0, sb = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	sa += adjointlambda(i, 0) * v(i, 0);
      // 	sb += adjointlambda(i, 1) * v(i, 0);
      // }
      // printf("t = %e, %e, s = %e, %e\n", ta, tb, sa, sb);
      
      dGvdp.resize(nparameters, 2);

      for (size_t k = 0; k < 2; k ++) {
	for (size_t j = 0; j < nparameters; j ++) {

	  
	  double a = 0.0;
	  for (size_t i = 0; i < size; i ++) {
	    
	    a += adjointlambda(i, k)*dCxv(i, j) * delta * gamma;
	    a += adjointlambda(size + i, k) * dCzv(i, j) * delta * gamma;
	    // a += adjointlambda(2*size + i, 0) * (omega*omega*dAxv(i, j) - dDxv(i, j)) * delta;
	    // a += adjointlambda(3*size + i, 0) * (omega*omega*dAzv(i, j) - dDzv(i, j)) * delta;

	    a += adjointlambda(2*size + i, k) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	    a += adjointlambda(3*size + i, k) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;
	    
	    a += sfundamental*adjointlambda(i, k)*dBxv(i, j) * delta * gamma * gamma;
	    a += sfundamental*adjointlambda(size + i, k)*dBzv(i, j) * delta * gamma * gamma;
	  }

	  dGvdp(j, k) = -a;
	}
      }
      
      eH = v(0, 0);
      eV = v(size, 0);
    }
    
    return k;
  }

  real solve_fundamental_gradient_v3(const Mesh<real, maxorder> &mesh,
				    size_t boundaryorder,
				    real omega,
				    Spec1DMatrix<real> &dgxdv,
				    Spec1DMatrix<real> &dgzdv,
				    Spec1DMatrix<real> &dkdp,
				    real &normA,
				    real &normB,
				    real &normC,
				    real &normD,
				    real &eH,
				    real &eV,
				    Spec1DMatrix<real> &dGvdp)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);

    adjointA.resize(4*size, 4*size);
    adjointB.resize(4*size, 1);
    adjointlambda.resize(4*size, 1);

    //
    // Precopy A matrix before it is destroyed by GEP
    //
    // for (size_t i = 0;  i < 4*size; i ++) {
    //   for (size_t j = 0; j < 4*size; j ++) {
    // 	adjointA(j, i) = As(i, j); // Note transpose
    //   }
    // }
    
    if (!SpecializedEigenProblem(As, Bs, work, eu, ev, lambda, Q, Z)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    real fundamental = 0.0;
    real sfundamental = 0.0;
    v.resize(4*size, 1);
    u.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);

	if (absk > fundamental) {
	  sfundamental = k;
	  fundamental = absk;

	  for (size_t j = 0; j < 4*size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    real k = 0.0;
    if (fundamental > 0.0) {
      //
      // Normalize u, v
      //
      double unorm = 0.0;
      double vnorm = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	unorm += u(i, 0) * u(i, 0);
	vnorm += v(i, 0) * v(i, 0);
      }
      unorm = sqrt(unorm);
      vnorm = sqrt(vnorm);
      for (size_t i = 0; i < 4*size; i ++) {
	u(i, 0) /= unorm;
	v(i, 0) /= vnorm;
      }
	
      size_t nparameters = postcomputegradient(mesh, boundaryorder, v);

      dkdp.resize(nparameters,1);
      dkdp.setZero();


      real norm = 0.0;
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      normD = 0.0;
      for (size_t j = 0; j < size; j ++) {
	norm -= u(j, 0) * (Bx(j, j) * gamma * gamma * delta) * v(j, 0) ;
	norm -= u(size + j, 0) * (Bz(j, j) * gamma * gamma * delta) * v(size + j, 0) ;
	norm -= u(2*size + j, 0) * v(2*size + j, 0);
	norm -= u(3*size + j, 0) * v(3*size + j, 0);

	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  v(j, 0) * Ax(j, j) * v(j, 0) +
	  v(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  v(j, 0) * Bx(j, j) * v(j, 0) +
	  v(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += v(j, 0) * cx + v(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += v(j, 0) * dx + v(size + j, 0) * dz;
      }

      k = sfundamental * gamma;

      for (size_t j = 0; j < nparameters; j ++) {

	double udAv = 0.0;
	double udBv = 0.0;
	double udCv = 0.0;
	double udDv = 0.0;
	
	for (size_t i = 0; i < size; i ++) {

	  dkdp(j, 0) += v(i, 0) * dCxv(i, j) * delta * gamma;
	  dkdp(j, 0) += v(i, 0) * sfundamental * dBxv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += v(size + i, 0) * dCzv(i, j) * delta * gamma;
	  dkdp(j, 0) += v(size + i, 0) * sfundamental * dBzv(i, j) * gamma * gamma * delta;

	  dkdp(j, 0) += v(2*size + i, 0) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  dkdp(j, 0) += v(3*size + i, 0) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;

	  udAv += u(i, 0)*dAxv(i, j) + u(size + i, 0)*dAzv(i, j);
	  udBv += u(i, 0)*dBxv(i, j) + u(size + i, 0)*dBzv(i, j);
	  udCv += u(i, 0)*dCxv(i, j) + u(size + i, 0)*dCzv(i, j);
	  udDv += u(i, 0)*dDxv(i, j) + u(size + i, 0)*dDzv(i, j);
	  
	}

	dkdp(j, 0) *= gamma/norm;

      }


      //
      // Adjoint matrix (initial set before call to GEP)
      //
      // for (size_t j = 0; j < size; j ++) {
      // 	adjointA(j, j) += Bx(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + size, j + size) += Bz(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + 2*size, j + 2*size) += sfundamental;
      // 	adjointA(j + 3*size, j + 3*size) += sfundamental;
      // }

      // sepdump("debugAs.txt", adjointA);
      
      // printf("Adjoint A\n");
      // for (size_t j = 0; j < 4*size; j ++) {
      // 	printf("%02d %16.9e\n", adjointA(j, 0));
      // }
      
      //
      // Projection of gradient to orthogonal direction.
      //
      eH = v(0, 0);
      eV = v(size, 0);
      adjointB.setZero();
      // printf("gp\n");
      real vdgx = 0.0;
      real vdgz = 0.0;
      for (size_t i = 0; i < size; i ++) {
	vdgx += v(i, 0) * dgxdv(i, 0)/eV;
	vdgz += v(size + i, 0) * dgzdv(i, 0) * (-eH/(eV*eV));
      }
      for (size_t i = 0; i < 4*size; i ++) {
	if (i < size) {
	  adjointB(i, 0) = dgxdv(i, 0)/eV;
	} else if (i < 2*size) {
	  adjointB(i, 0) = dgzdv(i - size, 0) * (-eH/(eV*eV));
	}
	
	adjointB(i, 0) -= v(i, 0) * (vdgx + vdgz);
      }

      // v.T adjointB == 0
      // real ta = 0.0, tb = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	ta += adjointB(i, 0) * v(i, 0);
      // 	tb += adjointB(i, 1) * v(i, 0);
      // }
      // SpecializedEigenProblemVerify(As, Q, Z, tAs, work);
      // SpecializedEigenProblemVerify(Bs, Q, Z, tBs, work);

      // printf("Check back\n");
      // double maxdiff = 0.0;
      // for (int i = 0; i < 4*size; i ++) {
      // 	for (int j = 0; j < 4*size; j ++) {

      // 	  double a = tAs(i, j) - sfundamental*tBs(i, j);
      // 	  double b = adjointA(j, i);

      // 	  double diff = fabs(a - b);
      // 	  if (diff > maxdiff) {
      // 	    maxdiff = diff;
      // 	  }
      // 	  if (diff > 1.0) {
      // 	    printf("  %2d, %2d : %16.9e\n", i, j, diff);
      // 	  }
      // 	}
      // }
      // printf("  max %16.9e\n", maxdiff );

      //
      // Solve for lambda_0 with forward substitution
      //
      adjointlambda.resize(4*size, 1);
      if (!SpecializedEigenProblemAdjoint<real>(As,
						Bs,
						Q,
						Z,
						sfundamental,
						adjointB,
						adjointlambda,
						work)) {
	ERROR("Failed to forward subsitute");
	return 0.0;
      }

      // printf("Lambda0\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	printf("%16.9e %16.9e\n", adjointlambda(i, 0), adjointlambda(i, 1));
      // }
      
      //
      // Check
      //
      // printf("Check a1 == b2 : a2 == b2 (using forward subs)\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real a = 0.0;
      // 	real b = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  a += adjointA(i, j) * adjointlambda(j, 0);
      // 	  b += adjointA(i, j) * adjointlambda(j, 1);
      // 	}
      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       a, adjointB(i, 0),
      // 	       b, adjointB(i, 1));
      // }

      //
      // Check
      //
      // printf("Check 2\n");
      // while (adjointA.rows() > IPIV_size) {
      // 	delete IPIV;
      // 	IPIV_size *= 2;
      // 	IPIV = new int[IPIV_size];
      // }
      // adjointAt.resize(adjointA.rows(), adjointA.cols());
      // adjointBt.resize(adjointB.rows(), adjointB.cols());
      
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  adjointAt(i, j) = adjointA(i, j);
      // 	}

      // 	for (size_t j = 0; j < adjointB.cols(); j ++) {
      // 	  adjointBt(i, j) = adjointB(i, j);
      // 	}
      // }
      
      // if (!GeneralSolve(adjointA, adjointB, IPIV)) {
      // 	FATAL("Failed to solve adjoint");
      // }
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	printf("%2d %16.9e %16.9e : %16.9e %16.9e\n",
      // 	       (int)i,
      // 	       adjointlambda(i, 0), adjointB(i, 0),
      // 	       adjointlambda(i, 1), adjointB(i, 1));
      // }

      //
      // Check 3
      //
      // printf("Check 3: a1 == b1 : a2 == b2 (using inverse)\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real a = 0.0;
      // 	real b = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  a += adjointAt(i, j) * adjointB(j, 0);
      // 	  b += adjointAt(i, j) * adjointB(j, 1);
      // 	}
      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       a, adjointBt(i, 0),
      // 	       b, adjointBt(i, 1));
      // }
      
      //
      // Ensure lambda0 orthogonal (Need to be orthogonal to u as non-symmetric)
      //
      real uA = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	uA += u(i, 0) * adjointlambda(i, 0);
      }
      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) -= u(i, 0)*uA;
      }

      // printf("Check 4: a1 == b1 : a2 == b2 projected\n");
      // uA = 0.0;
      // uB = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	uA += u(i, 0) * adjointB(i, 0);
      // 	uB += u(i, 0) * adjointB(i, 1);
      // }
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	adjointB(i, 0) -= u(i, 0)*uA;
      // 	adjointB(i, 1) -= u(i, 0)*uB;

      // 	printf("%16.9e %16.9e : %16.9e %16.9e\n",
      // 	       adjointlambda(i, 0), adjointB(i, 0),
      // 	       adjointlambda(i, 1), adjointB(i, 1));
      // }
      
      

      //
      // Adjoint lambda is now lambda_0, need to add correction for solution
      //
      //
      // w = vT B (B diagonal)
      //
      real wnorm = 0.0;
      adjointw.resize(4*size, 1);
      for (size_t i = 0; i < size; i ++) {

	adjointw(i, 0) = v(i, 0) * Bx(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i, 0)*adjointw(i, 0);

	adjointw(i + size, 0) = v(i + size, 0) * Bz(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i + size, 0)*adjointw(i + size, 0);

	adjointw(i + 2*size, 0) = v(i + 2*size, 0);
	wnorm += adjointw(i + 2*size, 0)*adjointw(i + 2*size, 0);

	adjointw(i + 3*size, 0) = v(i + 3*size, 0);
	wnorm += adjointw(i + 3*size, 0)*adjointw(i + 3*size, 0);

      }

      real dA_a = 0.0;
      real dB = 0.0;
      wnorm = sqrt(wnorm);
      
      for (size_t i = 0; i < 4*size; i ++) {
	adjointw(i, 0) /= wnorm;
	dA_a += adjointw(i, 0) * adjointlambda(i, 0);
	dB += adjointw(i, 0) * u(i, 0);
      }

      if (dB == 0.0) {
	FATAL("Adjoint w and u are orthogonal");
      }
      real adjointscale_a = -dA_a/dB;

      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) += adjointscale_a * u(i, 0);
      }

      // printf("A * lambda\n");
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	real c = 0.0;
      // 	for (size_t j = 0; j < 4*size; j ++) {
      // 	  c += adjointAt(i, j) * adjointlambda(j, 0);
      // 	}
      // 	printf("%02d %16.9e\n", i, c);
      // }
      
      // v.T adjointlambda == 0
      // real sa = 0.0, sb = 0.0;
      // for (size_t i = 0; i < 4*size; i ++) {
      // 	sa += adjointlambda(i, 0) * v(i, 0);
      // 	sb += adjointlambda(i, 1) * v(i, 0);
      // }
      // printf("t = %e, %e, s = %e, %e\n", ta, tb, sa, sb);
      
      dGvdp.resize(nparameters, 1);

      for (size_t j = 0; j < nparameters; j ++) {

	
	double a = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  
	  a += adjointlambda(i, k)*dCxv(i, j) * delta * gamma;
	  a += adjointlambda(size + i, k) * dCzv(i, j) * delta * gamma;
	  // a += adjointlambda(2*size + i, 0) * (omega*omega*dAxv(i, j) - dDxv(i, j)) * delta;
	  // a += adjointlambda(3*size + i, 0) * (omega*omega*dAzv(i, j) - dDzv(i, j)) * delta;
	  
	  a += adjointlambda(2*size + i, k) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  a += adjointlambda(3*size + i, k) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;
	  
	  a += sfundamental*adjointlambda(i, k)*dBxv(i, j) * delta * gamma * gamma;
	  a += sfundamental*adjointlambda(size + i, k)*dBzv(i, j) * delta * gamma * gamma;
	}
	
	dGvdp(j, k) = -a;
      }
      
      eH = v(0, 0);
      eV = v(size, 0);
    }
    
    return k;
  }

  real solve_fundamental_gradient_generic(const Mesh<real, maxorder> &mesh,
					  size_t boundaryorder,
					  real omega,
					  Spec1DMatrix<real> &dgxdv,
					  Spec1DMatrix<real> &dgzdv,
					  Spec1DMatrix<real> &dkdp,
					  Spec1DMatrix<real> &dUdp,
					  real &normA,
					  real &normB,
					  real &normC,
					  real &normD,
					  real &eH,
					  real &eV,
					  Spec1DMatrix<real> &dGvdp)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);

    adjointA.resize(4*size, 4*size);
    adjointB.resize(4*size, 3);
    adjointlambda.resize(4*size, 3);
    real galpha[3];
    
    //
    // Precopy A matrix before it is destroyed by GEP
    //
    // for (size_t i = 0;  i < 4*size; i ++) {
    //   for (size_t j = 0; j < 4*size; j ++) {
    // 	adjointA(j, i) = As(i, j); // Note transpose
    //   }
    // }
    
    if (!SpecializedEigenProblem(As, Bs, work, eu, ev, lambda, Q, Z)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    real fundamental = 0.0;
    real sfundamental = 0.0;
    v.resize(4*size, 1);
    u.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);

	if (k > 0.0 && absk > fundamental) {
	  sfundamental = k;
	  fundamental = absk;

	  for (size_t j = 0; j < 4*size; j ++) {
	    u(j, 0) = eu(j, i);
	    v(j, 0) = ev(j, i);
	  }
	}
      }
    }

    real k = 0.0;
    if (fundamental > 0.0) {
      //
      // Normalize u, v
      //
      double unorm = 0.0;
      double vnorm = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	unorm += u(i, 0) * u(i, 0);
	vnorm += v(i, 0) * v(i, 0);
      }
      unorm = sqrt(unorm);
      vnorm = sqrt(vnorm);
      for (size_t i = 0; i < 4*size; i ++) {
	u(i, 0) /= unorm;
	v(i, 0) /= vnorm;
      }
	
      size_t nparameters = postcomputegradient(mesh, boundaryorder, v);

      dkdp.resize(nparameters,1);
      dkdp.setZero();

      dUdp.resize(nparameters,1);
      dUdp.setZero();

      real norm = 0.0;
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      normD = 0.0;
      for (size_t j = 0; j < size; j ++) {
	norm -= u(j, 0) * (Bx(j, j) * gamma * gamma * delta) * v(j, 0) ;
	norm -= u(size + j, 0) * (Bz(j, j) * gamma * gamma * delta) * v(size + j, 0) ;
	norm -= u(2*size + j, 0) * v(2*size + j, 0);
	norm -= u(3*size + j, 0) * v(3*size + j, 0);

	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  v(j, 0) * Ax(j, j) * v(j, 0) +
	  v(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  v(j, 0) * Bx(j, j) * v(j, 0) +
	  v(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += v(j, 0) * cx + v(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += v(j, 0) * dx + v(size + j, 0) * dz;
      }

      k = sfundamental * gamma;

      //
      // Adjoint matrix (initial set before call to GEP)
      //
      // for (size_t j = 0; j < size; j ++) {
      // 	adjointA(j, j) += Bx(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + size, j + size) += Bz(j, j) * gamma * gamma * delta * sfundamental;
      // 	adjointA(j + 2*size, j + 2*size) += sfundamental;
      // 	adjointA(j + 3*size, j + 3*size) += sfundamental;
      // }

      // sepdump("debugAs.txt", adjointA);
      
      // printf("Adjoint A\n");
      // for (size_t j = 0; j < 4*size; j ++) {
      // 	printf("%02d %16.9e\n", adjointA(j, 0));
      // }
      
      //
      // Projection of gradient to orthogonal direction.
      //
      adjointB.setZero();
      // printf("gp\n");
      real vdgx = 0.0;
      real vdgz = 0.0;
      real vdgU = 0.0;
      for (size_t i = 0; i < size; i ++) {
	vdgx += v(i, 0) * dgxdv(i, 0);
	vdgz += v(size + i, 0) * dgzdv(i, 0);

	real av1 = Ax(i, i) * v(i, 0);
	real av2 = Az(i, i) * v(size + i, 0);
	real bv1 = Bx(i, i) * v(i, 0);
	real bv2 = Bz(i, i) * v(size + i, 0);
	real cvzz = 0.0;
	real cvxz = 0.0;
	real cvzx = 0.0;
	real cvxx = 0.0;
	  
	for (size_t j = 0; j < size; j ++) {
	  cvzz += Cz(j, i) * v(size + j, 0);
	  cvxz += Cx(i, j) * v(size + j, 0);
	  cvzx += Cz(i, j) * v(j, 0);
	  cvxx += Cx(j, i) * v(j, 0);
	}
	// printf("%2d %16.9e %16.9e %16.9e: %16.9e %16.9e %16.9e\n",
	//        i, cvzz, cvxz, cvzz + cvxz, cvzx, cvxx, cvzx + cvxx);
	
	adjointB(i, 2) =
	  (normA*(cvzz + cvxz) - 2.0*normC*av1 - 4.0*k*normB*av1 + 4.0*k*normA*bv1)/
	  (normA * normA * omega * 2.0);
	
	adjointB(size + i, 2) =
	  (normA*(cvzx + cvxx) - 2.0*normC*av2 - 4.0*k*normB*av2 + 4.0*k*normA*bv2)/
	  (normA * normA * omega * 2.0);

	vdgU +=
	  adjointB(i, 2) * v(i, 0) +
	  adjointB(size + i, 2) * v(size + i, 0);
      }

      //
      // Check dUdV against numerical approximation
      //
      // adjointw.resize(4*size, 1);
      // printf("norms: %16.9e %16.9e %16.9e\n", normA, normB, normC);
      // approximatedUdV(adjointw, 1.0e-4, k, omega);
      // printf("Check dUdV a ~= b\n");
      // for (size_t i = 0; i < 2*size; i ++) {
      // 	printf("%16.9e %16.9e\n", adjointB(i, 2), adjointw(i, 0));
      // }
      
      for (size_t i = 0; i < 4*size; i ++) {
	if (i < size) {
	  adjointB(i, 0) = dgxdv(i, 0);
	} else if (i < 2*size) {
	  adjointB(i, 1) = dgzdv(i - size, 0);
	}
	
	adjointB(i, 0) -= v(i, 0) * vdgx;
	adjointB(i, 1) -= v(i, 0) * vdgz;

	adjointB(i, 2) -= v(i, 0) * vdgU;

	// printf("%02d %16.9e\n", adjointB(i, 0));
      }

      galpha[0] = 0.0;
      galpha[1] = 0.0;
      galpha[2] = -(gamma * normB)/(omega * normA);

      //
      // Solve for lambda_0 with forward substitution
      //
      adjointlambda.resize(4*size, 3);
      if (!SpecializedEigenProblemAdjoint<real>(As,
						Bs,
						Q,
						Z,
						sfundamental,
						adjointB,
						adjointlambda,
						work)) {
	ERROR("Failed to forward subsitute");
	return 0.0;
      }
      
      //
      // Ensure lambda0 orthogonal (Need to be orthogonal to u as non-symmetric)
      //
      real uA = 0.0;
      real uB = 0.0;
      real uC = 0.0;
      for (size_t i = 0; i < 4*size; i ++) {
	uA += u(i, 0) * adjointlambda(i, 0);
	uB += u(i, 0) * adjointlambda(i, 1);
	uC += u(i, 0) * adjointlambda(i, 2);
      }
      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) -= u(i, 0)*uA;
	adjointlambda(i, 1) -= u(i, 0)*uB;
	adjointlambda(i, 2) -= u(i, 0)*uC;
      }

      //
      // Adjoint lambda is now lambda_0, need to add correction for solution
      //
      //
      // w = vT B (B diagonal)
      //
      real wnorm = 0.0;
      adjointw.resize(4*size, 1);
      for (size_t i = 0; i < size; i ++) {

	adjointw(i, 0) = v(i, 0) * Bx(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i, 0)*adjointw(i, 0);

	adjointw(i + size, 0) = v(i + size, 0) * Bz(i, i) * gamma*gamma*delta;
	wnorm += adjointw(i + size, 0)*adjointw(i + size, 0);

	adjointw(i + 2*size, 0) = v(i + 2*size, 0);
	wnorm += adjointw(i + 2*size, 0)*adjointw(i + 2*size, 0);

	adjointw(i + 3*size, 0) = v(i + 3*size, 0);
	wnorm += adjointw(i + 3*size, 0)*adjointw(i + 3*size, 0);

      }

      real dA_a = 0.0;
      real dA_b = 0.0;
      real dA_c = 0.0;
      real wdotu = 0.0;
      wnorm = sqrt(wnorm);
      
      for (size_t i = 0; i < 4*size; i ++) {
	adjointw(i, 0) /= wnorm;
	dA_a += adjointw(i, 0) * adjointlambda(i, 0);
	dA_b += adjointw(i, 0) * adjointlambda(i, 1);
	dA_c += adjointw(i, 0) * adjointlambda(i, 2);
	wdotu += adjointw(i, 0) * u(i, 0);
      }

      if (wdotu == 0.0) {
	FATAL("Adjoint w and u are orthogonal");
      }
      real adjointscale_a = (-dA_a - galpha[0]/wnorm)/wdotu;
      real adjointscale_b = (-dA_b - galpha[1]/wnorm)/wdotu;
      real adjointscale_c = (-dA_c - galpha[2]/wnorm)/wdotu;

      for (size_t i = 0; i < 4*size; i ++) {
	adjointlambda(i, 0) += adjointscale_a * u(i, 0);
	adjointlambda(i, 1) += adjointscale_b * u(i, 0);
	adjointlambda(i, 2) += adjointscale_c * u(i, 0);
      }

      //
      // Check -(B v)^T lambda = g_alpha
      //
      // for (int k = 0; k < 3; k ++) {
      // 	real bvTlambda = 0.0;
      // 	for (size_t i = 0; i < size; i ++) {
      // 	  bvTlambda += gamma*gamma*delta*Bx(i, i) * v(i, 0) * adjointlambda(i, k);
      // 	  bvTlambda += gamma*gamma*delta*Bz(i, i) * v(size + i, 0) * adjointlambda(size + i, k);
      // 	  bvTlambda += v(2*size + i, 0) * adjointlambda(2*size + i, k);
      // 	  bvTlambda += v(3*size + i, 0) * adjointlambda(3*size + i, k);
      // 	}
      // 	printf("  Check BvT Lambda %d: %16.9e %16.9e\n", k, bvTlambda, galpha[k]);
      // }
      
      dGvdp.resize(nparameters, 2);
      dUdp.resize(nparameters, 1);
      
      for (size_t l = 0; l < 3; l ++) {
	for (size_t j = 0; j < nparameters; j ++) {

	  
	  double a = 0.0;
	  double udAv = 0.0;
	  double udBv = 0.0;
	  double udCv = 0.0;
	  // double udDv = 0.0;
	  
	  for (size_t i = 0; i < size; i ++) {
	    
	    a += adjointlambda(i, l)*dCxv(i, j) * delta * gamma;
	    a += adjointlambda(size + i, l) * dCzv(i, j) * delta * gamma;

	    a += adjointlambda(2*size + i, l) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	    a += adjointlambda(3*size + i, l) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;
	    
	    a += sfundamental*adjointlambda(i, l)*dBxv(i, j) * delta * gamma * gamma;
	    a += sfundamental*adjointlambda(size + i, l)*dBzv(i, j) * delta * gamma * gamma;
	    
	    udAv += v(i, 0)*dAxv(i, j) + v(size + i, 0)*dAzv(i, j);
	    udBv += v(i, 0)*dBxv(i, j) + v(size + i, 0)*dBzv(i, j);
	    udCv += v(i, 0)*dCxv(i, j) + v(size + i, 0)*dCzv(i, j);
	    // udDv += v(i, 0)*dDxv(i, j) + v(size + i, 0)*dDzv(i, j);
	  }

	  if (l == 2) {
	    dUdp(j, 0) = (normA*(udCv + 2.0*k*udBv) -
			  udAv*(normC + 2.0*k*normB))/(2.0*omega*normA*normA) - a;	    
	  } else {
	    dGvdp(j, l) = -a;
	  }
	}
      }
      
      eH = v(0, 0);
      eV = v(size, 0);

      for (size_t j = 0; j < nparameters; j ++) {
	
	for (size_t i = 0; i < size; i ++) {
	  
	  dkdp(j, 0) += u(i, 0) * dCxv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(i, 0) * sfundamental * dBxv(i, j) * gamma * gamma * delta;
	  
	  dkdp(j, 0) += u(size + i, 0) * dCzv(i, j) * delta * gamma;
	  dkdp(j, 0) += u(size + i, 0) * sfundamental * dBzv(i, j) * gamma * gamma * delta;
	  
	  dkdp(j, 0) += u(2*size + i, 0) * (dDxv(i, j) - omega*omega*dAxv(i, j)) * delta;
	  dkdp(j, 0) += u(3*size + i, 0) * (dDzv(i, j) - omega*omega*dAzv(i, j)) * delta;
	  
	}
	
	dkdp(j, 0) *= gamma/norm;
	
      }
    }
    
    return k;
  }

  real solve_fundamental_scaled_vector(real omega,
				       const Mesh<real, maxorder> &mesh,
				       size_t boundaryorder,
				       MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude,
				       real &normA,
				       real &normB,
				       real &normC,
				       real &normD)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);

    if (!GEP(As, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    real fundamental = 0.0;
    real sfundamental = 0.0;

    u.resize(4*size, 1);
    v.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real k = lambda(i, 0)/lambda(i, 2);
	real absk = fabs(k);

	if (absk > fundamental) {
	  fundamental = absk;
	  sfundamental = k;

	  for (size_t j = 0; j < 4*size; j ++) {
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
      normD = 0.0;
      
      for (size_t j = 0; j < size; j ++) {
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  u(j, 0) * Ax(j, j) * v(j, 0) +
	  u(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  u(j, 0) * Bx(j, j) * v(j, 0) +
	  u(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += u(j, 0) * cx + u(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += u(j, 0) * dx + u(size + j, 0) * dz;
      }
    }

    fill_amplitude(mesh, v, boundaryorder, amplitude);

    return sfundamental * gamma;
  }
  
  bool solve_all_scaled_vector(real omega,
			       const Mesh<real, maxorder> &mesh,
			       size_t boundaryorder,
			       std::vector<real> &k,
			       std::vector<MeshAmplitude<real, maxorder, maxboundaryorder>> &amplitude,
			       real &normA,
			       real &normB,
			       real &normC,
			       real &normD)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);

    if (!GEP(As, Bs, work, eu, ev, lambda)) {
      ERROR("Failed to compute generalized eigen problem");
      return 0.0;
    }

    u.resize(4*size, 1);
    v.resize(4*size, 1);
    
    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {

	real ik = lambda(i, 0)/lambda(i, 2);

	k.push_back(ik * gamma);

	
	for (size_t j = 0; j < 4*size; j ++) {
	  u(j, 0) = eu(j, i);
	  v(j, 0) = ev(j, i);
	}

	size_t m = amplitude.size();
	amplitude.push_back(MeshAmplitude<real, maxorder, maxboundaryorder>());
	
	fill_amplitude(mesh, v, boundaryorder, amplitude[m]);
	
      }
    }
    
    if (k.size() > 0) {
      normA = 0.0;
      normB = 0.0;
      normC = 0.0;
      normD = 0.0;
      
      for (size_t j = 0; j < size; j ++) {
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	normA +=
	  u(j, 0) * Ax(j, j) * v(j, 0) +
	  u(size + j, 0) * Az(j, j) * v(size + j, 0);

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	normB += 
	  u(j, 0) * Bx(j, j) * v(j, 0) +
	  u(size + j, 0) * Bz(j, j) * v(size + j, 0);

	real cx = 0.0;
	real cz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  cx += Cx(j, i) * v(size + i, 0);
	  cz += Cz(j, i) * v(i, 0);
	}
	normC += u(j, 0) * cx + u(size + j, 0) * cz;

	real dx = 0.0;
	real dz = 0.0;
	for (size_t i = 0; i < size; i ++) {
	  dx += Dx(j, i) * v(i, 0);
	  dz += Dz(j, i) * v(size + i, 0);
	}
	normD += u(j, 0) * dx + u(size + j, 0) * dz;
      }
    }


    return k.size() > 0;
  }

  bool solve_all_scaled(real omega, std::vector<real> &k, real epsilon = 1.0e-9)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    real delta;
    real gamma = computeE_scaled(omega, delta);
    
    if (!GEP(As, Bs, work, eu, ev, lambda)) {
      FATAL("Failed to compute generalized eigen problem");
    }

    for (size_t i = 0; i < 4*size; i ++) {
      if (lambda(i, 1) == 0.0) {
	real absk = gamma * fabs(lambda(i, 0)/lambda(i, 2));
	bool unique = true;

	for (auto ek: k) {
	  if (fabs(ek - absk) < epsilon) {
	    unique = false;
	    break;
	  }
	}

	if (unique) {
	  k.push_back(absk);
	}
      }
    }

    return true;
  }

  void computeAx(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    Ax.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  Ax(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].rho *
	    Lobatto[mesh.cells[i].order]->weights[j];
	  
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  Ax(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].rho *
	    Lobatto[mesh.cells[i].order]->weights[j];
	  
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	Ax(offset + j, offset + j) +=
	  1.0/laguerrescalex *
	  mesh.boundary.rho *
	  Laguerre[boundaryorder]->weights[j];
      }
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	Ax(offset + j, offset + j) +=
	  mesh.cells[i].thickness/2.0 *
	  mesh.cells[i].nodes[j].rho *
	  Lobatto[mesh.cells[i].order]->weights[j];

      }
    }
  }

  void computeAxGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<real> &v)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    //
    // dAxv is non-zero for v(0 .. size - 1, 0) so we only calculate
    // that.
    dAxv.resize(size, nparameters);
    dAxv.setZero();

    if (mesh.cells.size() > 0) {

      for (size_t i = 0; i < (ncells - 1); i ++) {
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {

	    size_t l = mesh.cell_parameter_offsets[i] + k;
	  
	    dAxv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_RHO) *
	      v(offset + j, 0);

	  }
	  
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {

	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    
	    dAxv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_RHO) *
	      v(offset + j, 0);
	  }
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {

	size_t nparametersincell = mesh.boundary_jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.boundary_parameter_offset + k;
	  
	  dAxv(offset + j, l) +=
	    1.0/laguerrescalex *
	    Laguerre[boundaryorder]->weights[j] *
	    mesh.boundary_jacobian(k, MESHOFFSET_RHO) *
	    v(offset + j, 0);
	}
      }
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Row of dA/dp V, offset + j is col
	  dAxv(offset + j, l) +=
	    mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] *
	    mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_RHO) *
	    v(offset + j, 0);
	}

      }
    }
  }

  void computeAz(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    Az.setZero();

    if (mesh.cells.size() > 0) {
      
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	  Az(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].rho *
	    Lobatto[mesh.cells[i].order]->weights[j];
	  
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  Az(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].rho *
	    Lobatto[mesh.cells[i].order]->weights[j];
	  
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	Az(offset + j, offset + j) +=
	  1.0/laguerrescalez *
	  mesh.boundary.rho *
	  Laguerre[boundaryorder]->weights[j];
      }
      
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	Az(offset + j, offset + j) +=
	  mesh.cells[i].thickness/2.0 *
	  mesh.cells[i].nodes[j].rho *
	  Lobatto[mesh.cells[i].order]->weights[j];

      }
    }
  }

  void computeAzGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<real> &v)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    //
    // dAzv is non-zero for v(size .. 2*size - 1, 0) so we only calculate
    // that.
    dAzv.resize(size, nparameters);
    dAzv.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {

	    size_t l = mesh.cell_parameter_offsets[i] + k;
	  
	    dAzv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_RHO) *
	      v(size + offset + j, 0);

	  }
	  
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {

	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    
	    dAzv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_RHO) *
	      v(size + offset + j, 0);
	  }
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {

	size_t nparametersincell = mesh.boundary_jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.boundary_parameter_offset + k;
	  
	  dAzv(offset + j, l) +=
	    1.0/laguerrescalez *
	    Laguerre[boundaryorder]->weights[j] *
	    mesh.boundary_jacobian(k, MESHOFFSET_RHO) *
	    v(size + offset + j, 0);
	}
      }
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k; // Row of dA/dp V, offset + j is col
	  dAzv(offset + j, l) +=
	    mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] *
	    mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_RHO) *
	    v(size + offset + j, 0);
	}

      }
    }
  }

  void computeBx(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    Bx.setZero();

    if (mesh.cells.size() > 0) {

      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  Bx(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].A *
	    Lobatto[mesh.cells[i].order]->weights[j];
	  
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  Bx(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].A *
	    Lobatto[mesh.cells[i].order]->weights[j];
	  
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	Bx(offset + j, offset + j) +=
	  1.0/laguerrescalex *
	  mesh.boundary.A *
	  Laguerre[boundaryorder]->weights[j];
      }
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	Bx(offset + j, offset + j) +=
	  mesh.cells[i].thickness/2.0 *
	  mesh.cells[i].nodes[j].A *
	  Lobatto[mesh.cells[i].order]->weights[j];

      }
    }
      
  }

  void computeBxGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<real> &v)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    dBxv.resize(size, nparameters);
    dBxv.setZero();

    if (mesh.cells.size() > 0) {

      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    
	    dBxv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_A) *
	      v(offset + j, 0);
	  }
	    
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    dBxv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_A) *
	      v(offset + j, 0);
	  }
	  
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	size_t nparametersincell = mesh.boundary_jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.boundary_parameter_offset + k; // Row of dA/dp V, offset + j is col
	  dBxv(offset + j, l) +=
	    1.0/laguerrescalex *
	    Laguerre[boundaryorder]->weights[j] *
	    mesh.boundary_jacobian(k, MESHOFFSET_A) *
	    v(offset + j, 0);
	}
      }
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {
	
	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k;
	  dBxv(offset + j, l) +=
	    mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] *
	    mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_A) *
	    v(offset + j, 0);
	}
      }
    }
  }

  void computeBz(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    Bz.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  Bz(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].L *
	    Lobatto[mesh.cells[i].order]->weights[j];

	}

	offset += mesh.cells[i].order;
      }
    }
    
    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  Bz(offset + j, offset + j) +=
	    mesh.cells[i].thickness/2.0 *
	    mesh.cells[i].nodes[j].L *
	    Lobatto[mesh.cells[i].order]->weights[j];
	  
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	Bz(offset + j, offset + j) +=
	  1.0/laguerrescalez *
	  mesh.boundary.L *
	  Laguerre[boundaryorder]->weights[j];
      }
      
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	Bz(offset + j, offset + j) +=
	  mesh.cells[i].thickness/2.0 *
	  mesh.cells[i].nodes[j].L *
	  Lobatto[mesh.cells[i].order]->weights[j];

      }
    }
  }

  void computeBzGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<real> &v)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    dBzv.resize(size, nparameters);
    dBzv.setZero();

    if (mesh.cells.size() > 0) {

      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    
	    dBzv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_L) *
	      v(size + offset + j, 0);
	  }
	    
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    dBzv(offset + j, l) +=
	      mesh.cells[i].thickness/2.0 *
	      Lobatto[mesh.cells[i].order]->weights[j] *
	      mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_L) *
	      v(size + offset + j, 0);
	  }
	  
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	size_t nparametersincell = mesh.boundary_jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.boundary_parameter_offset + k; // Row of dA/dp V, offset + j is col
	  dBzv(offset + j, l) +=
	    1.0/laguerrescalex *
	    Laguerre[boundaryorder]->weights[j] *
	    mesh.boundary_jacobian(k, MESHOFFSET_L) *
	    v(size + offset + j, 0);
	}
      }
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {
	
	size_t nparametersincell = mesh.cells[i].jacobian.rows();
	for (size_t k = 0; k < nparametersincell; k ++) {
	  size_t l = mesh.cell_parameter_offsets[i] + k;
	  dBzv(offset + j, l) +=
	    mesh.cells[i].thickness/2.0 *
	    Lobatto[mesh.cells[i].order]->weights[j] *
	    mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_L) *
	    v(size + offset + j, 0);
	}
      }
    }
  }

  void computeCx(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    Cx.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  

	    Cx(offset + m, offset + j) +=
	      (Lobatto[mesh.cells[i].order]->weights[m] * 
	       mesh.cells[i].nodes[m].F *
	       Lobatto[mesh.cells[i].order]->derivative_weights[j][m] -
	       
	       Lobatto[mesh.cells[i].order]->weights[j] *
	       mesh.cells[i].nodes[j].L *
	       Lobatto[mesh.cells[i].order]->derivative_weights[m][j]);
	  }
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      // Do the last Lobatto cell
      if (mesh.cells.size() > 0) {
	size_t i = mesh.cells.size() - 1;
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	    
	    Cx(offset + m, offset + j) +=
	      (Lobatto[mesh.cells[i].order]->weights[m] * 
	       mesh.cells[i].nodes[m].F *
	       Lobatto[mesh.cells[i].order]->derivative_weights[j][m] -
	       
	       Lobatto[mesh.cells[i].order]->weights[j] *
	       mesh.cells[i].nodes[j].L *
	       Lobatto[mesh.cells[i].order]->derivative_weights[m][j]);
	  }

	}
	
	offset += mesh.cells[i].order;
      }
      
      
      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t j = 0; j <= boundaryorder; j ++) {

	  
	  Cx(offset + m, offset + j) +=
	    (Laguerre[boundaryorder]->weights[m] * 
	     mesh.boundary.F *
	     Laguerre[boundaryorder]->derivative_weights[j][m] -

	     Laguerre[boundaryorder]->weights[j] *
	     mesh.boundary.L *
	     Laguerre[boundaryorder]->derivative_weights[m][j]);
	}
      }
      
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	  
	  Cx(offset + m, offset + j) +=
	    (Lobatto[mesh.cells[i].order]->weights[m] * 
	     mesh.cells[i].nodes[m].F *
	     Lobatto[mesh.cells[i].order]->derivative_weights[j][m] -

	     Lobatto[mesh.cells[i].order]->weights[j] *
	     mesh.cells[i].nodes[j].L *
	     Lobatto[mesh.cells[i].order]->derivative_weights[m][j]);
	}
      }
    }
  }

  void computeCxGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<double> &v)
  {
    dCxv.resize(size, nparameters);
    dCxv.setZero();
    
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    dCxv.setZero();
    
    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;

	      dCxv(offset + m, l) += 
		(Lobatto[mesh.cells[i].order]->weights[m] * 
		 mesh.cells[i].jacobian(k, m*6 + MESHOFFSET_F) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[j][m] -
	       
		 Lobatto[mesh.cells[i].order]->weights[j] *
		 mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_L) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[m][j]) * v(size + offset + j, 0);
	    }
	  }
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      // Do the last Lobatto cell
      if (mesh.cells.size() > 0) {
	size_t i = mesh.cells.size() - 1;
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  
	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;

	      dCxv(offset + m, l) += 
		(Lobatto[mesh.cells[i].order]->weights[m] * 
		 mesh.cells[i].jacobian(k, m*6 + MESHOFFSET_F) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[j][m] -
	       
		 Lobatto[mesh.cells[i].order]->weights[j] *
		 mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_L) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[m][j]) * v(size + offset + j, 0);
	    }
	  }

	}
	
	offset += mesh.cells[i].order;
      }
      
      
      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t j = 0; j <= boundaryorder; j ++) {
	  
	  size_t nparametersincell = mesh.boundary_jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.boundary_parameter_offset + k;
	    dCxv(offset + m, l) +=
	      (Laguerre[boundaryorder]->weights[m] * 
	       mesh.boundary_jacobian(k, MESHOFFSET_F) *
	       Laguerre[boundaryorder]->derivative_weights[j][m] -
	       
	       Laguerre[boundaryorder]->weights[j] *
	       mesh.boundary_jacobian(k, MESHOFFSET_L) *
	       Laguerre[boundaryorder]->derivative_weights[m][j]) *
	      v(size + offset + j, 0);
	  }	      
	}
      }
      
    } else {
      
      size_t i = mesh.cells.size() - 1;
      
      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	for (size_t j = 0; j < mesh.cells[i].order; j ++) {
	  
	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    
	    dCxv(offset + m, l) += 
	      (Lobatto[mesh.cells[i].order]->weights[m] * 
	       mesh.cells[i].jacobian(k, m*6 + MESHOFFSET_F) *
	       Lobatto[mesh.cells[i].order]->derivative_weights[j][m] -
	       
	       Lobatto[mesh.cells[i].order]->weights[j] *
	       mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_L) *
	       Lobatto[mesh.cells[i].order]->derivative_weights[m][j]) * v(size + offset + j, 0);
	  }
	}
      }
    }
  }

  void computeCz(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    Cz.setZero();

    if (mesh.cells.size() > 0) {

      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	    Cz(offset + m, offset + j) +=
	      (Lobatto[mesh.cells[i].order]->weights[j] * 
	       mesh.cells[i].nodes[j].F *
	       Lobatto[mesh.cells[i].order]->derivative_weights[m][j] -
	       
	       Lobatto[mesh.cells[i].order]->weights[m] *
	       mesh.cells[i].nodes[m].L *
	       Lobatto[mesh.cells[i].order]->derivative_weights[j][m]);
	  }
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last Lobatto cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	    
	    Cz(offset + m, offset + j) +=
	      (Lobatto[mesh.cells[i].order]->weights[j] * 
	       mesh.cells[i].nodes[j].F *
	       Lobatto[mesh.cells[i].order]->derivative_weights[m][j] -
	       
	       Lobatto[mesh.cells[i].order]->weights[m] *
	       mesh.cells[i].nodes[m].L *
	       Lobatto[mesh.cells[i].order]->derivative_weights[j][m]);
	  }
	  
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t j = 0; j <= boundaryorder; j ++) {

	  
	  Cz(offset + m, offset + j) +=
	    (Laguerre[boundaryorder]->weights[j] * 
	     mesh.boundary.F *
	     Laguerre[boundaryorder]->derivative_weights[m][j] -

	     Laguerre[boundaryorder]->weights[m] *
	     mesh.boundary.L *
	     Laguerre[boundaryorder]->derivative_weights[j][m]);
	}
      }
      
    } else {
      size_t i = mesh.cells.size() - 1;
      
      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	for (size_t j = 0; j < mesh.cells[i].order; j ++) {

	  
	  Cz(offset + m, offset + j) +=
	    (Lobatto[mesh.cells[i].order]->weights[j] * 
	     mesh.cells[i].nodes[j].F *
	     Lobatto[mesh.cells[i].order]->derivative_weights[m][j] -

	     Lobatto[mesh.cells[i].order]->weights[m] *
	     mesh.cells[i].nodes[m].L *
	     Lobatto[mesh.cells[i].order]->derivative_weights[j][m]);
	}
      }
    }
  }

  void computeCzGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<double> &v)
  {
    dCzv.resize(size, nparameters);
    dCzv.setZero();
    
    size_t ncells = mesh.cells.size();
    size_t offset = 0;

    dCzv.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  

	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;

	      dCzv(offset + m, l) += 
		(Lobatto[mesh.cells[i].order]->weights[j] * 
		 mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_F) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[m][j] -
	       
		 Lobatto[mesh.cells[i].order]->weights[m] *
		 mesh.cells[i].jacobian(k, m*6 + MESHOFFSET_L) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[j][m]) * v(offset + j, 0);
	    }
	  }
	}

	offset += mesh.cells[i].order;
      }
    }

    if (mesh.boundary.rho > 0.0) {

      //
      // Laguerre
      //

      // Do the last Lobatto cell
      if (mesh.cells.size() > 0) {
	size_t i = mesh.cells.size() - 1;
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	    
	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;

	      dCzv(offset + m, l) += 
		(Lobatto[mesh.cells[i].order]->weights[j] * 
		 mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_F) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[m][j] -
	       
		 Lobatto[mesh.cells[i].order]->weights[m] *
		 mesh.cells[i].jacobian(k, m*6 + MESHOFFSET_L) *
		 Lobatto[mesh.cells[i].order]->derivative_weights[j][m]) * v(offset + j, 0);
	    }
	  }

	}
	
	offset += mesh.cells[i].order;
      }
      
      
      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t j = 0; j <= boundaryorder; j ++) {
	  
	  size_t nparametersincell = mesh.boundary_jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.boundary_parameter_offset + k;
	    dCzv(offset + m, l) +=
	      (Laguerre[boundaryorder]->weights[j] * 
	       mesh.boundary_jacobian(k, MESHOFFSET_F) *
	       Laguerre[boundaryorder]->derivative_weights[m][j] -
	       
	       Laguerre[boundaryorder]->weights[m] *
	       mesh.boundary_jacobian(k, MESHOFFSET_L) *
	       Laguerre[boundaryorder]->derivative_weights[j][m]) *
	      v(offset + j, 0);
	  }	      
	}
      }
      
    } else {
      
      size_t i = mesh.cells.size() - 1;
      
      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	for (size_t j = 0; j < mesh.cells[i].order; j ++) {
	  
	  size_t nparametersincell = mesh.cells[i].jacobian.rows();
	  for (size_t k = 0; k < nparametersincell; k ++) {
	    size_t l = mesh.cell_parameter_offsets[i] + k;
	    
	    dCzv(offset + m, l) += 
	      (Lobatto[mesh.cells[i].order]->weights[j] * 
	       mesh.cells[i].jacobian(k, j*6 + MESHOFFSET_F) *
	       Lobatto[mesh.cells[i].order]->derivative_weights[m][j] -
	       
	       Lobatto[mesh.cells[i].order]->weights[m] *
	       mesh.cells[i].jacobian(k, m*6 + MESHOFFSET_L) *
	       Lobatto[mesh.cells[i].order]->derivative_weights[j][m]) * v(offset + j, 0);
	  }
	}
      }
    }
  }

  void computeDx(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    Dx.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	    for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	      
	      Dx(offset + m, offset + n) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].nodes[j].L * 
		Lobatto[mesh.cells[i].order]->weights[j] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		Lobatto[mesh.cells[i].order]->derivative_weights[m][j];
	      
	    }
	  }
	}

	offset += mesh.cells[i].order;
      }
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	    for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
		
		Dx(offset + m, offset + n) += 
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
      
      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t n = 0; n <= boundaryorder; n ++) {
	  for (size_t j = 0; j <= boundaryorder; j ++) {
	      
	    Dx(offset + m, offset + n) += 
	      laguerrescalex *
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

      // Note subtle difference in use of < instead of <= as fixed boundary is zero displacement so
      // we remove it from the loop.
      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	for (size_t n = 0; n < mesh.cells[i].order; n ++) { 
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	    Dx(offset + m, offset + n) += 
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

  void computeDxGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<real> &v)
  {
    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    dDxv.resize(size, nparameters);
    dDxv.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	    for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	      size_t nparametersincell = mesh.cells[i].jacobian.rows();
	      
	      for (size_t k = 0; k < nparametersincell; k ++) {
		size_t l = mesh.cell_parameter_offsets[i] + k;
		dDxv(offset + m, l) += 
		  2.0/mesh.cells[i].thickness *
		  mesh.cells[i].jacobian(k, 6*j + MESHOFFSET_L) * 
		  Lobatto[mesh.cells[i].order]->weights[j] * 
		  Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		  Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		  v(offset + n, 0);
	      }
	    }
	  }
	}

	offset += mesh.cells[i].order;
      }
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	    for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	      size_t nparametersincell = mesh.cells[i].jacobian.rows();
	      
	      for (size_t k = 0; k < nparametersincell; k ++) {
		size_t l = mesh.cell_parameter_offsets[i] + k;
		dDxv(offset + m, l) += 
		  2.0/mesh.cells[i].thickness *
		  mesh.cells[i].jacobian(k, 6*j + MESHOFFSET_L) * 
		  Lobatto[mesh.cells[i].order]->weights[j] * 
		  Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		  Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		  v(offset + n, 0);
	      }
	    }
	  }
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t n = 0; n <= boundaryorder; n ++) {
	  for (size_t j = 0; j <= boundaryorder; j ++) {

	    size_t nparametersincell = mesh.boundary_jacobian.rows();

	    for (size_t k = 0; k < nparametersincell; k ++) {

	      size_t l = mesh.boundary_parameter_offset + k;
	      
	      dDxv(offset + m, l) += 
		laguerrescalex *
		mesh.boundary_jacobian(k, MESHOFFSET_L) * 
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

      // Note subtle difference in use of < instead of <= as fixed boundary is zero displacement so
      // we remove it from the loop.
      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	for (size_t n = 0; n < mesh.cells[i].order; n ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	    
	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;
	      dDxv(offset + m, l) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].jacobian(k, 6*j + MESHOFFSET_L) * 
		Lobatto[mesh.cells[i].order]->weights[j] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		v(offset + n, 0);
	    }
	  }
	}
      }
    }
  }

  void computeDz(const Mesh<real, maxorder> &mesh, size_t boundaryorder)
  {
    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    Dz.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  for (size_t k = 0; k <= mesh.cells[i].order; k ++) {
	    for (size_t l = 0; l <= mesh.cells[i].order; l ++) {
	      
	      Dz(offset + j, offset + k) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].nodes[l].C * 
		Lobatto[mesh.cells[i].order]->weights[l] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[j][l] *
		Lobatto[mesh.cells[i].order]->derivative_weights[k][l];
	      
	    }
	  }
	}

	offset += mesh.cells[i].order;
      }
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	  for (size_t k = 0; k <= mesh.cells[i].order; k ++) {
	    for (size_t l = 0; l <= mesh.cells[i].order; l ++) {
	      
	      Dz(offset + j, offset + k) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].nodes[l].C * 
		Lobatto[mesh.cells[i].order]->weights[l] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[j][l] *
		Lobatto[mesh.cells[i].order]->derivative_weights[k][l];
		
	    }
	  }
	}
	
	offset += mesh.cells[i].order;
      }

      // Do the Laguerre cell
      for (size_t j = 0; j <= boundaryorder; j ++) {
	for (size_t k = 0; k <= boundaryorder; k ++) {
	  for (size_t l = 0; l <= boundaryorder; l ++) {
	      
	    Dz(offset + j, offset + k) += 
	      laguerrescalez *
	      mesh.boundary.C * 
	      Laguerre[boundaryorder]->weights[l] * 
	      Laguerre[boundaryorder]->derivative_weights[j][l] *
	      Laguerre[boundaryorder]->derivative_weights[k][l];

	  }
	}
      }

    } else {
      //
      // Fixed
      //
      size_t i = ncells - 1;

      for (size_t j = 0; j < mesh.cells[i].order; j ++) {
	// Note subtle difference in use of < instead of <= as fixed boundary is zero displacement so
	// we remove it from the loop.
	for (size_t k = 0; k < mesh.cells[i].order; k ++) { 
	  for (size_t l = 0; l <= mesh.cells[i].order; l ++) {

	    Dz(offset + j, offset + k) += 
	      2.0/mesh.cells[i].thickness *
	      mesh.cells[i].nodes[l].C * 
	      Lobatto[mesh.cells[i].order]->weights[l] * 
	      Lobatto[mesh.cells[i].order]->derivative_weights[j][l] *
	      Lobatto[mesh.cells[i].order]->derivative_weights[k][l];
	  }
	}
      }
    }
  }

  void computeDzGradient(size_t nparameters,
			 const Mesh<real, maxorder> &mesh,
			 size_t boundaryorder,
			 Spec1DMatrix<real> &v)
  {
    size_t ncells = mesh.cells.size();

    size_t offset = 0;
    dDzv.resize(size, nparameters);
    dDzv.setZero();

    if (mesh.cells.size() > 0) {
      for (size_t i = 0; i < (ncells - 1); i ++) {
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	    for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	      size_t nparametersincell = mesh.cells[i].jacobian.rows();
	      
	      for (size_t k = 0; k < nparametersincell; k ++) {
		size_t l = mesh.cell_parameter_offsets[i] + k;
		dDzv(offset + m, l) += 
		  2.0/mesh.cells[i].thickness *
		  mesh.cells[i].jacobian(k, 6*j + MESHOFFSET_C) * 
		  Lobatto[mesh.cells[i].order]->weights[j] * 
		  Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		  Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		  v(size + offset + n, 0);
	      }
	    }
	  }
	}

	offset += mesh.cells[i].order;
      }
    }
    
    //
    // Last cell is different depending on Boundary condition
    //
    if (mesh.boundary.rho > 0.0) {
      //
      // Laguerre
      //

      if (mesh.cells.size() > 0) {
	
	// Do the last cell
	size_t i = mesh.cells.size() - 1;
	
	for (size_t m = 0; m <= mesh.cells[i].order; m ++) {
	  for (size_t n = 0; n <= mesh.cells[i].order; n ++) {
	    for (size_t j = 0; j <= mesh.cells[i].order; j ++) {

	      size_t nparametersincell = mesh.cells[i].jacobian.rows();
	      
	      for (size_t k = 0; k < nparametersincell; k ++) {
		size_t l = mesh.cell_parameter_offsets[i] + k;
		dDzv(offset + m, l) += 
		  2.0/mesh.cells[i].thickness *
		  mesh.cells[i].jacobian(k, 6*j + MESHOFFSET_C) * 
		  Lobatto[mesh.cells[i].order]->weights[j] * 
		  Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		  Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		  v(size + offset + n, 0);
	      }
	    }
	  }
	}
	
	offset += mesh.cells[i].order;
      }
      
      // Do the Laguerre cell
      for (size_t m = 0; m <= boundaryorder; m ++) {
	for (size_t n = 0; n <= boundaryorder; n ++) {
	  for (size_t j = 0; j <= boundaryorder; j ++) {

	    size_t nparametersincell = mesh.boundary_jacobian.rows();

	    for (size_t k = 0; k < nparametersincell; k ++) {

	      size_t l = mesh.boundary_parameter_offset + k;
	      
	      dDzv(offset + m, l) += 
		laguerrescalex *
		mesh.boundary_jacobian(k, MESHOFFSET_C) * 
		Laguerre[boundaryorder]->weights[j] * 
		Laguerre[boundaryorder]->derivative_weights[m][j] *
		Laguerre[boundaryorder]->derivative_weights[n][j] *
		v(size + offset + n, 0);
	    }
	  }
	}
      }

    } else {
      //
      // Fixed
      //
      size_t i = ncells - 1;

      // Note subtle difference in use of < instead of <= as fixed boundary is zero displacement so
      // we remove it from the loop.
      for (size_t m = 0; m < mesh.cells[i].order; m ++) {
	for (size_t n = 0; n < mesh.cells[i].order; n ++) {
	  for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	    
	    size_t nparametersincell = mesh.cells[i].jacobian.rows();
	    
	    for (size_t k = 0; k < nparametersincell; k ++) {
	      size_t l = mesh.cell_parameter_offsets[i] + k;
	      dDzv(offset + m, l) += 
		2.0/mesh.cells[i].thickness *
		mesh.cells[i].jacobian(k, 6*j + MESHOFFSET_C) * 
		Lobatto[mesh.cells[i].order]->weights[j] * 
		Lobatto[mesh.cells[i].order]->derivative_weights[m][j] *
		Lobatto[mesh.cells[i].order]->derivative_weights[n][j] *
		v(size + offset + n, 0);
	    }
	  }
	}
      }
    }
  }

  void computeE(real omega)
  {
    //
    // E is the first companion form for solving the quadratic eigen value problem
    //
    // E = | I  0 |^-1 | 0                I |
    //     | 0 -B |    | D - omega^2 A    C |
    //
    // with
    // A = | Ax  0 |
    //     |  0 Az |
    // B = | Bx  0 |
    //     |  0 Bz |
    // C = |  0 Cx |
    //     | Cz  0 |
    // D = | Dx  0 |
    //     |  0 Dz |
    //

    if (size == 0) {
      FATAL("Uninitialized");
    }
    
    real o2 = omega * omega;

    E.setZero();

    //
    // Top left quadrant all zero
    //
    
    //
    // Top right quadrant identity
    //
    for (size_t j = 0; j < size; j ++) {
      E(j, 2*size + j) = 1.0;
      E(size + j, 3*size + j) = 1.0;
    }
    
    //
    // Bottom left quadrant (omega^2*A - D)/Bx
    //
    for (size_t j = 0; j < size; j ++) {

      for (size_t i = 0; i < size; i ++) {

	E(2*size + j, i) = (o2*Ax(j, i) - Dx(j, i))/Bx(j, j);
	E(3*size + j, size + i) = (o2*Az(j, i) - Dz(j, i))/Bz(j, j);

      }
    }

    //
    // Bottom right quadrant C
    //
    for (size_t j = 0; j < size; j ++) {

      for (size_t i = 0; i < size; i ++) {

	E(2*size + j, 3*size + i) = -Cx(j, i)/Bx(j, j);
	E(3*size + j, 2*size + i) = -Cz(j, i)/Bz(j, j);

      }
    }
    
  }

  real computeE_scaled(real omega, real &_delta)
  {
    real o2 = omega * omega;

    A0x.setZero();
    A0z.setZero();

    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {

	A0x(j, i) = Dx(j, i) - Ax(j, i)*o2;
	A0z(j, i) = Dz(j, i) - Az(j, i)*o2;

      }
    }

    real gamma = 1.0;
    real delta = 1.0;

    if (!disable_scale) {
      double a = A0x.norm();
      double b = A0z.norm();
      double A0norm = sqrt(a*a + b*b);
      
      a = Cx.norm();
      b = Cz.norm();
      double A1norm = sqrt(a*a + b*b);
      
      a = Bx.norm();
      b = Bz.norm();
      double A2norm = sqrt(a*a + b*b);
      
      gamma = sqrt(A0norm/A2norm);
      delta = 2.0/(A0norm + gamma * A1norm);
    }
    
    As.setZero();
    Bs.setZero();
    
    for (size_t j = 0; j < size; j ++) {
      //
      // -I
      //
      As(j, j + 2*size) = -1.0;
      As(j + size, j + 3*size) = -1.0;

      Bs(j + 2*size, j + 2*size) = -1.0;
      Bs(j + 3*size, j + 3*size) = -1.0;
      
      //
      // A2
      //
      Bs(j, j) = -gamma*gamma*delta * Bx(j, j);
      Bs(j + size, j + size) = -gamma*gamma*delta * Bz(j, j);
      
      for (size_t i = 0; i < size; i ++) {

	//
	// A1
	//
	As(j, i + size) = gamma * delta * Cx(j, i);
	As(j + size, i) = gamma * delta * Cz(j, i);

	//
	// A0
	//
	As(j + 2*size, i) = delta * A0x(j, i);
	As(j + 3*size, i + size) = delta * A0z(j, i);

      }
    }

    _delta = delta;
    return gamma;
  }

  void computeEx(real omega)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    Ex.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        Ex(j, i) = (o2 * Ax(j, i) - Dx(j, i))/Bx(j, j);
      }
    }
  }

  void computeEz(real omega)
  {
    if (size == 0) {
      FATAL("Unconfigured");
    }

    Ez.setZero();

    real o2 = omega*omega;
    for (size_t j = 0; j < size; j ++) {
      for (size_t i = 0; i < size; i ++) {
        Ez(j, i) = (o2 * Az(j, i) - Dz(j, i))/Bz(j, j);
      }
    }
  }

  real computeBlockSquareMatrixProductSum(Spec1DMatrix<real> &Mx,
					  Spec1DMatrix<real> &Mz,
					  Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < size; i ++) {
      real sx = 0.0;
      real sz = 0.0;
      
      for (size_t j = 0; j < size; j ++) {
	sx += v(i, 0) * Mx(j, i);
	sz += v(size + i, 0) * Mz(j, i);
      }

      s += sx*v(i, 0) + sz*v(size + i, 0);
    }

    return s;
  }

  real computeReverseBlockSquareMatrixProductSum(Spec1DMatrix<real> &Mx,
						 Spec1DMatrix<real> &Mz,
						 Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < size; i ++) {
      real sx = 0.0;
      real sz = 0.0;
      
      for (size_t j = 0; j < size; j ++) {
	sx += v(i, 0) * Mx(j, i);
	sz += v(size + i, 0) * Mz(j, i);
      }

      s += sx*v(size + i, 0) + sz*v(i, 0);
    }

    return s;
  }

  real computeDiagonalBlockSquareMatrixProductSum(Spec1DMatrix<real> &Mx,
						  Spec1DMatrix<real> &Mz,
						  Spec1DMatrix<double> &v)
  {
    real s = 0.0;
    
    for (size_t i = 0; i < size; i ++) {
      s +=
	v(i, 0) * Mx(i, i) * v(i, 0) +
	v(size + i, 0) * Mz(i, i) * v(size + i, 0);
    }

    return s;
  }

  real computeLinearSquareMatrixProductSum(Spec1DMatrix<double> &u,
					   Spec1DMatrix<real> &a,
					   Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < 4*size; i ++) {
      real n = 0.0;
      for (size_t j = 0; j < 4*size; j ++) {
	n += u(j, 0) * a(j, i);
      }
      
      s += n * v(i, 0);
    }

    return s;
  }

  real computeLinearSquareDiagonalMatrixProductSum(Spec1DMatrix<double> &u,
						   real &gamma,
						   real &delta,
						   Spec1DMatrix<real> &Bx,
						   Spec1DMatrix<real> &Bz,
						   Spec1DMatrix<double> &v)
  {
    real s = 0.0;

    for (size_t i = 0; i < size; i ++) {
      s -= u(i, 0) * gamma*gamma*delta*Bx(i, i) * v(i, 0);
      s -= u(i + size, 0) * gamma*gamma*delta*Bz(i, i) * v(i + size, 0);
      s -= u(i + 2*size, 0) * v(i + 2*size, 0);
      s -= u(i + 3*size, 0) * v(i + 3*size, 0);
    }

    return s;
  }

  real computeLinearSquareMatrixProductSumAsymmetric(Spec1DMatrix<double> &u,
						     real &gamma,
						     real &delta,
						     Spec1DMatrix<real> &Dx,
						     Spec1DMatrix<real> &Dz,
						     Spec1DMatrix<double> &v)
  {
    real s1 = 0.0;
    real s2 = 0.0;

    for (size_t i = 0; i < size; i ++) {

      real n1 = 0.0;
      real n2 = 0.0;
      for (size_t j = 0; i < size; i ++) {

	n1 += delta * Dx(j, i) * u(2*size + i, 0);
	n2 += delta * Dz(j, i) * u(3*size + i, 0);

      }

      s1 += n1 * v(i, 0);
      s2 += n2 * v(size + i, 0);
    }

    return s1 + s2;
  }

  real computeGroupSquareMatrixProduct(Spec1DMatrix<double> &u,
				       Spec1DMatrix<real> &X,
				       Spec1DMatrix<real> &Z,
				       Spec1DMatrix<double> &v)
  {
    real s = 0.0;
    for (size_t i = 0; i < size; i ++) {
      real n1 = 0.0;
      real n2 = 0.0;
      for (size_t j = 0; j < size; j ++) {
	n1 += u(j, 0) * X(j, i);
	n2 += u(size + j, 0) * Z(j, i);
      }

      s += n1*v(i, 0) + n2*v(size + i, 0);
    }
    return s;
  }
				       
  real computeGroupAsymmetricSquareMatrixProduct(Spec1DMatrix<double> &u,
						 Spec1DMatrix<real> &X,
						 Spec1DMatrix<real> &Z,
						 Spec1DMatrix<double> &v)
  {
    real s = 0.0;
    for (size_t i = 0; i < size; i ++) {
      real n1 = 0.0;
      real n2 = 0.0;
      for (size_t j = 0; j < size; j ++) {
	n1 += u(j + size, 0) * Z(j, i);
	n2 += u(j, 0) * X(j, i);
      }

      s += n1*v(i, 0) + n2*v(size + i, 0);
    }
    return s;
  }

  real computeGroupSquareDiagonalMatrixProduct(Spec1DMatrix<double> &u,
					       Spec1DMatrix<real> &X,
					       Spec1DMatrix<real> &Z,
					       Spec1DMatrix<double> &v)
  {
    real s = 0.0;
    for (size_t i = 0; i < size; i ++) {
      s += u(i, 0) * X(i, i) * v(i, 0) +
	u(size + i, 0) * Z(i, i) * v(size + i, 0);
    }
    return s;
  }

  void printMatrix(FILE *fp, const Spec1DMatrix<real> &m) const
  {
    if (size > 0) {
      for (size_t j = 0; j < size; j ++) {
	for (size_t i = 0; i < size; i ++) {
	  
	  fprintf(fp, "%10.3g ", (double)m(j, i));
	}
	fprintf(fp, "\n");
      }
    }
  }
  
  void printAx(FILE *fp) const
  {
    printMatrix(fp, Ax);
  }
  void printBx(FILE *fp) const
  {
    printMatrix(fp, Bx);
  }
  void printCx(FILE *fp) const
  {
    printMatrix(fp, Cx);
  }
  void printDx(FILE *fp) const
  {
    printMatrix(fp, Dx);
  }
  
  void printAz(FILE *fp) const
  {
    printMatrix(fp, Az);
  }
  void printBz(FILE *fp) const
  {
    printMatrix(fp, Bz);
  }
  void printCz(FILE *fp) const
  {
    printMatrix(fp, Cz);
  }
  void printDz(FILE *fp) const
  {
    printMatrix(fp, Dz);
  }

  void printE(FILE *fp) const
  {
    if (size > 0) {
      for (size_t j = 0; j < (4*size); j ++) {
	for (size_t i = 0; i < (4*size); i ++) {
	  
	  fprintf(fp, "%.15g ", E(j, i));
	}
	fprintf(fp, "\n");
      }
    }
  }
  
  void printELLQx(FILE *fp) const
  {
    if (size > 0) {
      for (size_t j = 0; j < size; j ++) {
	for (size_t i = 0; i < size; i ++) {
	  
	  fprintf(fp, "%10.3g ", E(j + 2*size, i));
	}
	fprintf(fp, "\n");
      }
    }
  }

  void printELLQz(FILE *fp) const
  {
    if (size > 0) {
      for (size_t j = 0; j < size; j ++) {
	for (size_t i = 0; i < size; i ++) {
	  
	  fprintf(fp, "%10.3g ", E(j + 3*size, size + i));
	}
	fprintf(fp, "\n");
      }
    }
  }

  void printELRQx(FILE *fp) const
  {
    if (size > 0) {
      for (size_t j = 0; j < size; j ++) {
	for (size_t i = 0; i < size; i ++) {
	  
	  fprintf(fp, "%10.3g ", E(j + 2*size, 3*size + i));
	}
	fprintf(fp, "\n");
      }
    }
  }

  void printELRQz(FILE *fp) const
  {
    if (size > 0) {
      for (size_t j = 0; j < size; j ++) {
	for (size_t i = 0; i < size; i ++) {
	  
	  fprintf(fp, "%10.3g ", E(j + 3*size, 2*size + i));
	}
	fprintf(fp, "\n");
      }
    }
  }

  void fill_amplitude(const Mesh<real, maxorder> &mesh,
		      const Spec1DMatrix<real> &v,
		      size_t boundaryorder,
		      MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude)
  {
    int offset = 0;

    amplitude.cells.clear();
    
    for (size_t i = 0; i < mesh.cells.size() - 1; i ++) {

      CellAmplitude<real, maxorder> ca;

      ca.order = mesh.cells[i].order;
      ca.thickness = mesh.cells[i].thickness;
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	ca.nodes[j].V = 0.0;
	
	ca.nodes[j].U = v(offset + j, 0);
	ca.nodes[j].W = v(size + offset + j, 0);
      }

      amplitude.cells.push_back(ca);

      offset += mesh.cells[i].order;
    }

    if (mesh.boundary.rho == 0.0) {

      //
      // Last cell
      //
      size_t i = mesh.cells.size() - 1;
      
      CellAmplitude<real, maxorder> ca;

      ca.order = mesh.cells[i].order;
      ca.thickness = mesh.cells[i].thickness;
      for (size_t j = 0; j < mesh.cells[i].order; j ++) {
	ca.nodes[j].V = 0.0;
	ca.nodes[j].U = v(offset + j, 0);
	ca.nodes[j].W = v(size + offset + j, 0);
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
      for (size_t j = 0; j <= mesh.cells[i].order; j ++) {
	ca.nodes[j].V = 0.0;
	ca.nodes[j].U = v(offset + j, 0);
	ca.nodes[j].W = v(size + offset + j, 0);
      }

      amplitude.cells.push_back(ca);

      offset += mesh.cells[i].order;

      //
      // Laguerre halfspace
      //
      amplitude.boundary.thickness = 1.0;
      amplitude.boundary.order = boundaryorder;
      for (size_t j = 0; j <= boundaryorder; j ++) {
	amplitude.boundary.nodes[j].V = 0.0;
	amplitude.boundary.nodes[j].U = v(offset + j, 0);
	amplitude.boundary.nodes[j].W = v(size + offset + j, 0);
      }
    }
  }
  
  void interpolate_eigenvector(const MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude,
			       real z,
			       real &u, real &w)
  {
    real dz = z;
    for (auto &c : amplitude.cells) {

      if (dz < c.thickness) {

	real xi = 2.0*dz/c.thickness - 1.0;

	if (xi < -1.0 || xi > 1.0) {
	  FATAL("xi out of range: %f", xi);
	}
	
	u = 0.0;
	w = 0.0;

	for (size_t i = 0; i <= c.order; i ++) {

	  real s = Lobatto[c.order]->cardinal[i].value(xi);

	  u += s*c.nodes[i].U;
	  w += s*c.nodes[i].W;

	}

	return;

      } else {

	dz -= c.thickness;

      }
    }

    // Traversed cells, point in boundary/halfspace

    if (amplitude.boundary.thickness == 0.0) {
      //
      // Fixed boundary -> return 0
      //
      u = 0.0;
      w = 0.0;

    } else {
      //
      // Laguerre halfspace: dz = distance into halfspace, xi = dz*scale
      //
      real xiu = dz * laguerrescalex;
      real xiw = dz * laguerrescalez;

      u = 0.0;
      w = 0.0;
      
      for (size_t i = 0; i <= amplitude.boundary.order; i ++) {
	u += Laguerre[amplitude.boundary.order]->cardinal[i].value(xiu)*amplitude.boundary.nodes[i].U;
	w += Laguerre[amplitude.boundary.order]->cardinal[i].value(xiw)*amplitude.boundary.nodes[i].W;
      }
      
    }
  }

  void interpolate_eigenvector_derived(const Mesh<real, maxorder> &mesh,
				       const MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude,
				       real z,
				       real k,
				       real &r2, real &r4)
  {
    
    MeshParameter<real> p = mesh.interpolate(z);

    real dz = z;
    for (auto &c : amplitude.cells) {

      if (dz < c.thickness) {

	real xi = 2.0*dz/c.thickness - 1.0;

	if (xi < -1.0 || xi > 1.0) {
	  FATAL("xi out of range: %f", xi);
	}
	
	r2 = 0.0;
	r4 = 0.0;

	for (size_t i = 0; i <= c.order; i ++) {

	  real s = Lobatto[c.order]->cardinal[i].value(xi);
	  real ds = Lobatto[c.order]->cardinal[i].dydx(xi) * 2.0/c.thickness;

	  r2 += p.L*(ds*c.nodes[i].U - k*s*c.nodes[i].W);
	  r4 += p.F*k*s*c.nodes[i].U + p.C*ds*c.nodes[i].W;

	}

	return;

      } else {

	dz -= c.thickness;

      }
    }

    // Traversed cells, point in boundary/halfspace

    if (amplitude.boundary.thickness == 0.0) {
      //
      // Fixed boundary -> return 0
      //
      r2 = 0.0;
      r4 = 0.0;

    } else {
      //
      // Laguerre halfspace: dz = distance into halfspace, xi = dz*scale
      //
      real xiu = dz * laguerrescalex;
      real xiw = dz * laguerrescalez;

      r2 = 0.0;
      r4 = 0.0;
      
      for (size_t i = 0; i <= amplitude.boundary.order; i ++) {

	real s = Laguerre[amplitude.boundary.order]->cardinal[i].value(xiu);
	real ds = Laguerre[amplitude.boundary.order]->cardinal[i].dydx(xiu) * laguerrescalex;
	  
	r2 += p.L*(ds*amplitude.boundary.nodes[i].U - k*s*amplitude.boundary.nodes[i].W);
	
	s = Laguerre[amplitude.boundary.order]->cardinal[i].value(xiw);
	ds = Laguerre[amplitude.boundary.order]->cardinal[i].dydx(xiw) * laguerrescalez;
	  
	r4 += p.F*k*s*amplitude.boundary.nodes[i].U + p.C*ds*amplitude.boundary.nodes[i].W;
      }
      
    }
  }

  void get_eigenvector_nodes(const Mesh<real, maxorder> &mesh,
			     const MeshAmplitude<real, maxorder, maxboundaryorder> &amplitude,
			     real k,
			     std::vector<real> &z12,
			     std::vector<real> &r1,
			     std::vector<real> &r2,
			     std::vector<real> &z34,
			     std::vector<real> &r3,
			     std::vector<real> &r4)
  {
    real zoffset;
    size_t ci;
    
    zoffset = 0.0;
    ci = 0;
    for (auto &c : amplitude.cells) {

      for (size_t i = 0; i <= c.order; i ++) {
	
	real xi = Lobatto[c.order]->nodes[i];
	
	z12.push_back(zoffset + (xi + 1.0)/2.0 * c.thickness);
	z34.push_back(zoffset + (xi + 1.0)/2.0 * c.thickness);

	r1.push_back(c.nodes[i].U);
	r3.push_back(c.nodes[i].W);
	
	real r2s = -k * c.nodes[i].W * mesh.cells[ci].nodes[i].L;
	real r4s = mesh.cells[ci].nodes[i].F * k * c.nodes[i].U;

	//
	// Derivative components
	//
	for (size_t j = 0; j <= c.order; j ++) {
	  real ds = Lobatto[c.order]->cardinal[j].dydx(xi) * 2.0/c.thickness;
	  r2s += mesh.cells[ci].nodes[i].L * c.nodes[j].U * ds;
	  r4s += mesh.cells[ci].nodes[i].C * c.nodes[j].W * ds;
	}
	r2.push_back(r2s);
	r4.push_back(r4s);
	
      }

      zoffset += c.thickness;
      ci ++;
    }

    // Traversed cells, point in boundary/halfspace

    if (amplitude.boundary.thickness > 0.0) {
      //
      // Laguerre halfspace: dz = distance into halfspace, xi = dz*scale
      //
      
      for (size_t i = 0; i <= amplitude.boundary.order; i ++) {

	real xi = Laguerre[amplitude.boundary.order]->nodes[i];
	
	z12.push_back(zoffset + xi/laguerrescalex);
	r1.push_back(amplitude.boundary.nodes[i].U);
	
	real r2s = -k * amplitude.boundary.nodes[i].W * mesh.boundary.L;
	for (size_t j = 0; j <= amplitude.boundary.order; j ++) {
	  real ds = Laguerre[amplitude.boundary.order]->cardinal[j].dydx(xi) * laguerrescalex;
	  r2s += mesh.boundary.L * amplitude.boundary.nodes[j].U * ds;
	}
	r2.push_back(r2s);

	z34.push_back(zoffset + xi/laguerrescalez);
	r3.push_back(amplitude.boundary.nodes[i].W);

	//
	// Derivative components
	//
	real r4s = mesh.boundary.F * k * amplitude.boundary.nodes[i].U;
	for (size_t j = 0; j <= amplitude.boundary.order; j ++) {
	  real ds = Laguerre[amplitude.boundary.order]->cardinal[j].dydx(xi) * laguerrescalez;
	  r4s += mesh.cells[ci].nodes[i].C * amplitude.boundary.nodes[j].W * ds;
	}
	r4.push_back(r4s);
      }
      
    }
  }
  
  void approximatedUdV(Spec1DMatrix<double> &dUdv, double epsilon, double k, double omega)
  {
    for (int i = 0; i < 2*size; i ++) {

      double normA1 = 0.0;
      double normB1 = 0.0;
      double normC1 = 0.0;
      double normD1 = 0.0;

      //
      // U + dU
      //
      for (size_t j = 0; j < size; j ++) {
	if (j == i) {
	  normA1 +=
	    (v(j, 0) + epsilon/2.0) * Ax(j, j) * (v(j, 0) + epsilon/2.0) +
	    v(size + j, 0) * Az(j, j) * v(size + j, 0);
	} else if (size + j == i) {
	  normA1 +=
	    v(j, 0) * Ax(j, j) * v(j, 0) +
	    (v(size + j, 0) + epsilon/2.0) * Az(j, j) * (v(size + j, 0) + epsilon/2.0);
	} else {
	  normA1 +=
	    v(j, 0) * Ax(j, j) * v(j, 0) +
	    v(size + j, 0) * Az(j, j) * v(size + j, 0);
	}

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	if (j == i) {
	  normB1 += 
	    (v(j, 0) + epsilon/2.0) * Bx(j, j) * (v(j, 0) + epsilon/2.0) +
	    v(size + j, 0) * Bz(j, j) * v(size + j, 0);
	} else if (size + j == i) {
	  normB1 += 
	    v(j, 0) * Bx(j, j) * v(j, 0) +
	    (v(size + j, 0) + epsilon/2.0) * Bz(j, j) * (v(size + j, 0) + epsilon/2.0);
	} else {
	  normB1 += 
	    v(j, 0) * Bx(j, j) * v(j, 0) +
	    v(size + j, 0) * Bz(j, j) * v(size + j, 0);
	}

	real cx = 0.0;
	real cz = 0.0;
	for (size_t l = 0; l < size; l ++) {

	  if (size + l == i) {
	    cx += Cx(j, l) * (v(size + l, 0) + epsilon/2.0);
	  } else {
	    cx += Cx(j, l) * v(size + l, 0);
	  }

	  if (l == i) {
	    cz += Cz(j, l) * (v(l, 0) + epsilon/2.0);
	  } else {
	    cz += Cz(j, l) * v(l, 0);
	  }
	}

	if (j == i) {
	  normC1 += (v(j, 0) + epsilon/2.0) * cx;
	} else {
	  normC1 += v(j, 0) * cx;
	}

	if (size + j == i) {
	  normC1 += (v(size + j, 0) + epsilon/2.0) * cz;
	} else {
	  normC1 += v(size + j, 0) * cz;
	}

      }

      // printf("dnorm: %16.9e %16.9e %16.9e :", normA1, normB1, normC1);
      double U1 = (2.0*normB1*k + normC1)/(2.0*omega*normA1);

      double normA2 = 0.0;
      double normB2 = 0.0;
      double normC2 = 0.0;

      //
      // U - dU
      //
      for (size_t j = 0; j < size; j ++) {
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)
	// normA is: A diagonal
	//   u(2n .. 3n - 1) * Ax * v(0 .. n - 1)
	//   u(3n .. 4n - 1) * Az * v(n .. 2n - 1)

	if (j == i) {
	  normA2 +=
	    (v(j, 0) - epsilon/2.0) * Ax(j, j) * (v(j, 0) - epsilon/2.0) +
	    v(size + j, 0) * Az(j, j) * v(size + j, 0);
	} else if (size + j == i) {
	  normA2 +=
	    v(j, 0) * Ax(j, j) * v(j, 0) +
	    (v(size + j, 0) - epsilon/2.0) * Az(j, j) * (v(size + j, 0) - epsilon/2.0);
	} else {
	  normA2 +=
	    v(j, 0) * Ax(j, j) * v(j, 0) +
	    v(size + j, 0) * Az(j, j) * v(size + j, 0);
	}

	// normB is: B diagonal
	//   u(0 .. n - 1) * Bx * v(0 .. n - 1)
	//   u(n .. 2n - 1) * Bz * v(n .. 2n - 1)
	if (j == i) {
	  normB2 += 
	    (v(j, 0) - epsilon/2.0) * Bx(j, j) * (v(j, 0) - epsilon/2.0) +
	    v(size + j, 0) * Bz(j, j) * v(size + j, 0);
	} else if (size + j == i) {
	  normB2 += 
	    v(j, 0) * Bx(j, j) * v(j, 0) +
	    (v(size + j, 0) - epsilon/2.0) * Bz(j, j) * (v(size + j, 0) - epsilon/2.0);
	} else {
	  normB2 += 
	    v(j, 0) * Bx(j, j) * v(j, 0) +
	    v(size + j, 0) * Bz(j, j) * v(size + j, 0);
	}

	real cx = 0.0;
	real cz = 0.0;
	for (size_t l = 0; l < size; l ++) {

	  if (size + l == i) {
	    cx += Cx(j, l) * (v(size + l, 0) - epsilon/2.0);
	  } else {
	    cx += Cx(j, l) * v(size + l, 0);
	  }

	  if (l == i) {
	    cz += Cz(j, l) * (v(l, 0) - epsilon/2.0);
	  } else {
	    cz += Cz(j, l) * v(l, 0);
	  }
	}

	if (j == i) {
	  normC2 += (v(j, 0) - epsilon/2.0) * cx;
	} else {
	  normC2 += v(j, 0) * cx;
	}

	if (size + j == i) {
	  normC2 += (v(size + j, 0) - epsilon/2.0) * cz;
	} else {
	  normC2 += v(size + j, 0) * cz;
	}

      }

      double U2 = (2.0*normB2*k + normC2)/(2.0*omega*normA2);
      // printf(" %16.9e %16.9e %16.9e (%16.9e %16.9e)\n", normA2, normB2, normC2, U1, U2);
      // if (i < size) {
      // 	printf("dB/dV %16.9e %16.9e\n", (normB1 - normB2)/epsilon, 2.0*v(i, 0)*Bx(i, i));
      // } else {
      // 	printf("dB/dV %16.9e %16.9e\n", (normB1 - normB2)/epsilon, 2.0*v(i, 0)*Bz(i - size, i - size));
      // }

      if (i < size) {
	double vzcz = 0.0;
	double cxvz = 0.0;
	for (int j = 0; j < size; j ++) {
	  vzcz += v(size + j, 0) * Cz(j, i);
	  cxvz += Cx(i, j) * v(size + j, 0);
	} 
	printf("dC/dV %02d %16.9e %16.9e (%16.9e %16.9e)\n", i, (normC1 - normC2)/epsilon, vzcz + cxvz,
	       vzcz, cxvz);
      } else {
	double vxcx = 0.0;
	double czvx = 0.0;
	for (int j = 0; j < size; j ++) {
	  vxcx += v(j, 0) * Cx(j, i - size);
	  czvx += Cz(i - size, j) * v(j, 0);
	} 
	printf("dC/dV %02d %16.9e %16.9e (%16.9e %16.9e)\n", i, (normC1 - normC2)/epsilon, vxcx + czvx,
	       vxcx, czvx);
      }	
      
	       
      dUdv(i, 0) = (U1 - U2)/epsilon;

      
    }
  }
  
  size_t size;

  std::array<LobattoQuadrature<double, maxorder>*, maxorder + 1> Lobatto;
  std::array<LaguerreQuadrature<double, maxboundaryorder>*, maxboundaryorder + 1> Laguerre;

  real laguerrescalex;
  real laguerrescalez;

  Spec1DMatrix<real> A0x, A0z;

  Spec1DMatrix<real> Ax, Az;
  Spec1DMatrix<real> Bx, Bz;
  Spec1DMatrix<real> Cx, Cz;
  Spec1DMatrix<real> Dx, Dz;

  Spec1DMatrix<real> E;

  Spec1DMatrix<real> Ex;
  Spec1DMatrix<real> Ez;

  Spec1DMatrix<real> As;
  Spec1DMatrix<real> Bs;

  Spec1DMatrix<double> dAs;
  Spec1DMatrix<double> dBs;

  Spec1DMatrix<real> eu, ev;
  Spec1DMatrix<real> lambda;
  Spec1DMatrix<real> work;

  Spec1DMatrix<double> deu, dev;
  Spec1DMatrix<double> dlambda;
  Spec1DMatrix<double> dwork;

  Spec1DMatrix<real> u, v;

  Spec1DMatrix<real> dAxv, dAzv;
  Spec1DMatrix<real> dBxv, dBzv;
  Spec1DMatrix<real> dCxv, dCzv;
  Spec1DMatrix<real> dDxv, dDzv;
  
  Spec1DMatrix<real> Q;
  Spec1DMatrix<real> Z;

  Spec1DMatrix<real> adjointA;
  Spec1DMatrix<real> adjointAt;
  Spec1DMatrix<real> adjointB;
  Spec1DMatrix<real> adjointBt;
  Spec1DMatrix<real> adjointlambda;
  Spec1DMatrix<real> adjointw;

  Spec1DMatrix<real> tAs;
  Spec1DMatrix<real> tBs;
  
  int *IPIV;
  int IPIV_size;

  bool disable_scale;
};

#endif // rayleighmatrices_hpp


