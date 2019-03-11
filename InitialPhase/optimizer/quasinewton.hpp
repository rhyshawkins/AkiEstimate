#pragma once
#ifndef quasinewton_hpp
#define quasinewton_hpp

#include "common.hpp"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

//
// Quasi-Newton step (see eqn 3.59 of Tarantola:2005:A)
//
// m_{n+1} = m_n - \mu (G^T C_D^{-1} G + C_M^{-1})^{-1} (G^T C_D^{-1} (d_n - d_{obs}) + C_M^{-1}(m_n - m_0))
//

class QuasiNewton : public LeastSquaresIterator {
public:

  QuasiNewton() :
    A(nullptr),
    GTCdinv(nullptr),
    y(nullptr),
    mnp1(nullptr),
    perm(nullptr)
  {
  }
  
  virtual bool ComputeStep(double epsilon,
			   Spec1DMatrix<double> &C_d,
			   Spec1DMatrix<double> &C_m,
			   Spec1DMatrix<double> &residuals,
			   Spec1DMatrix<double> &G,
			   Spec1DMatrix<double> &dLdp,
			   Spec1DMatrix<int> &model_mask,
			   Spec1DMatrix<double> &current_model,
			   Spec1DMatrix<double> &prior_model,
			   Spec1DMatrix<double> &proposed_model)
  {
    int Nd = residuals.rows();
    int Nm = current_model.rows();
    
    allocate(Nd, Nm);

    //
    // Set G^T C_d^-1 (assuming Cd is diagonal, ie just a vector)
    //
    gsl_matrix_set_zero(GTCdinv);
    for (int i = 0; i < Nm; i ++) {
      for (int j = 0; j < Nd; j ++) {

	double s = G(j, i)/C_d(j, 0);
	gsl_matrix_set(GTCdinv, i, j, s);

      }
    }
	  
    
    //
    // Set A (we assume C_d and C_m are diagonal covariance matrices)
    //
    gsl_matrix_set_zero(A);
    for (int i = 0; i < Nm; i ++) {
      for (int j = 0; j < Nm; j ++) {

	double s = 0.0;
	if (i == j) {
	  s = 1.0/C_m(i, 0);
	}
	  
	for (int k = 0; k < Nd; k ++) {
	  s += gsl_matrix_get(GTCdinv, i, k) * G(k, j);
	}

	gsl_matrix_set(A, i, j, s);
      }

      // printf("\n");
    }

    //
    // Set y vector (A m_current - mu [G^T Cd^-1 residuals + Cm^-1 (m_current - m_0)])
    //
    gsl_vector_set_zero(y);

    for (int i = 0; i < Nm; i ++) {

      double s = -epsilon * (current_model(i, 0) - prior_model(i, 0))/C_m(i, 0);

      for (int j = 0; j < Nm; j ++) {

	s += gsl_matrix_get(A, i, j) * current_model(j, 0);
      }

      for (int j = 0; j < Nd; j ++) {
      	s -= epsilon * (gsl_matrix_get(GTCdinv, i, j) * residuals(j, 0));
      }

      gsl_vector_set(y, i, s);
    }

    //
    // Solve A m_{n + 1} = y
    //
    int signum = 0;
    if (gsl_linalg_LU_decomp(A, perm, &signum) < 0) {
      fprintf(stderr, "error: failed to LU decomp matrix\n");
      return false;
    }

    if (gsl_linalg_LU_solve(A, perm, y, mnp1) < 0) {
      fprintf(stderr, "error: failed to LU solve new position\n");
      return false;
    }

    //
    // Copy new model back
    //
    for (int i = 0; i < Nm; i ++) {
      proposed_model(i, 0) = gsl_vector_get(mnp1, i);
      // printf("%2d %16.9e\n", i, proposed_model(i, 0));
    }

    printf(" %16.9e -> %16.9e: %16.9e -> %16.9e: %16.9e -> %16.9e: %16.9e -> %16.9e\n\n",
	   current_model(0, 0), proposed_model(0, 0),
	   current_model(6, 0), proposed_model(6, 0),
	   current_model(12, 0), proposed_model(12, 0),
	   current_model(18, 0), proposed_model(18, 0));
    
    return true;
  }

  virtual bool ComputeStepJoint(double epsilon,
				Spec1DMatrix<double> &C_d_love,
				Spec1DMatrix<double> &C_d_rayleigh,
				Spec1DMatrix<double> &C_m,
				Spec1DMatrix<double> &residuals_love,
				Spec1DMatrix<double> &residuals_rayleigh,
				Spec1DMatrix<double> &G_love,
				Spec1DMatrix<double> &G_rayleigh,
				Spec1DMatrix<double> &dLdp,
				Spec1DMatrix<int> &model_mask,
				Spec1DMatrix<double> &current_model,
				Spec1DMatrix<double> &prior_model,
				Spec1DMatrix<double> &proposed_model)
  {
    int Nd_love = residuals_love.rows();
    int Nd_rayleigh = residuals_rayleigh.rows();
    int Nd = Nd_love + Nd_rayleigh;
    int Nm = current_model.rows();
    
    allocate(Nd, Nm);

    //
    // Set G^T C_d^-1 (assuming Cd is diagonal, ie just a vector)
    //
    gsl_matrix_set_zero(GTCdinv);
    for (int i = 0; i < Nm; i ++) {
      for (int j = 0; j < Nd_love; j ++) {

	double s = G_love(j, i)/C_d_love(j, 0);
	gsl_matrix_set(GTCdinv, i, j, s);

      }

      for (int j = 0; j < Nd_rayleigh; j ++) {

	double s = G_rayleigh(j, i)/C_d_rayleigh(j, 0);
	gsl_matrix_set(GTCdinv, i, Nd_love + j, s);

      }
    }
	  
    
    //
    // Set A (we assume C_d and C_m are diagonal covariance matrices)
    //
    gsl_matrix_set_zero(A);
    for (int i = 0; i < Nm; i ++) {
      for (int j = 0; j < Nm; j ++) {

	double s = 0.0;
	if (i == j) {
	  s = 1.0/C_m(i, 0);
	}
	  
	for (int k = 0; k < Nd_love; k ++) {
	  s += gsl_matrix_get(GTCdinv, i, k) * G_love(k, j);
	}

	for (int k = 0; k < Nd_rayleigh; k ++) {
	  s += gsl_matrix_get(GTCdinv, i, Nd_love + k) * G_rayleigh(k, j);
	}

	gsl_matrix_set(A, i, j, s);
      }

      // printf("\n");
    }

    //
    // Set y vector (A m_current - mu [G^T Cd^-1 residuals + Cm^-1 (m_current - m_0)])
    //
    gsl_vector_set_zero(y);

    for (int i = 0; i < Nm; i ++) {

      double s = -epsilon * (current_model(i, 0) - prior_model(i, 0))/C_m(i, 0);

      for (int j = 0; j < Nm; j ++) {

	s += gsl_matrix_get(A, i, j) * current_model(j, 0);
      }

      for (int j = 0; j < Nd_love; j ++) {
      	s -= epsilon * (gsl_matrix_get(GTCdinv, i, j) * residuals_love(j, 0));
      }

      for (int j = 0; j < Nd_rayleigh; j ++) {
      	s -= epsilon * (gsl_matrix_get(GTCdinv, i, Nd_love + j) * residuals_rayleigh(j, 0));
      }

      gsl_vector_set(y, i, s);
    }

    //
    // Solve A m_{n + 1} = y
    //
    int signum = 0;
    if (gsl_linalg_LU_decomp(A, perm, &signum) < 0) {
      fprintf(stderr, "error: failed to LU decomp matrix\n");
      return false;
    }

    if (gsl_linalg_LU_solve(A, perm, y, mnp1) < 0) {
      fprintf(stderr, "error: failed to LU solve new position\n");
      return false;
    }

    //
    // Copy new model back
    //
    for (int i = 0; i < Nm; i ++) {
      proposed_model(i, 0) = gsl_vector_get(mnp1, i);
      // printf("%2d %16.9e\n", i, proposed_model(i, 0));
    }

    printf(" %16.9e -> %16.9e: %16.9e -> %16.9e: %16.9e -> %16.9e: %16.9e -> %16.9e\n\n",
	   current_model(0, 0), proposed_model(0, 0),
	   current_model(6, 0), proposed_model(6, 0),
	   current_model(12, 0), proposed_model(12, 0),
	   current_model(18, 0), proposed_model(18, 0));
    
    return true;
  }

  
  void allocate(int Nd, int Nm)
  {
    if (A == nullptr || (int)A->size1 != Nm) {
      if (A != nullptr) {
	gsl_matrix_free(A);
      }
      A = gsl_matrix_alloc(Nm, Nm);
    }

    if (GTCdinv == nullptr || (int)GTCdinv->size1 != Nm || (int)GTCdinv->size2 != Nd) {
      if (GTCdinv != nullptr) {
	gsl_matrix_free(GTCdinv);
      }
      GTCdinv = gsl_matrix_alloc(Nm, Nd);
    }

    if (perm == nullptr || (int)perm->size != Nm) {
      if (perm != nullptr) {
	gsl_permutation_free(perm);
      }
      perm = gsl_permutation_alloc(Nm);
    }

    if (y == nullptr || (int)y->size != Nm) {
      if (y != nullptr) {
	gsl_vector_free(y);
      }
      y = gsl_vector_alloc(Nm);
    }

    if (mnp1 == nullptr || (int)mnp1->size != Nm) {
      if (mnp1 != nullptr) {
	gsl_vector_free(mnp1);
      }
      mnp1 = gsl_vector_alloc(Nm);
    }
    
  }


  gsl_matrix *A;
  gsl_matrix *GTCdinv;

  gsl_vector *y;
  gsl_vector *mnp1;
  
  gsl_permutation *perm;
};

#endif // quasinewton_hpp
  
