#ifndef common_hpp
#define common_hpp

#include "spec1d/anisotropicrhovsxivpvs.hpp"
#include "spec1d/anisotropicrhovsxivpvshalfspace.hpp"
#include "spec1d/empiricalmodel.hpp"
#include "spec1d/model.hpp"
#include "spec1d/modelloader.hpp"
#include "spec1d/lovematrices.hpp"
#include "spec1d/rayleighmatrices.hpp"
#include "spec1d/lobattoprojection.hpp"

constexpr int MAXORDER = 20;
constexpr int BOUNDARYORDER = MAXORDER;

constexpr double EPSILON_MIN = 1.0e-9;

typedef AnisotropicRhoVsXiVpVs<double> cell_parameter_t;
typedef AnisotropicRhoVsXiVpVsHalfspace<double, MAXORDER> boundary_t;
typedef Cell<double, cell_parameter_t, MAXORDER> cell_t;
typedef Model<double, cell_parameter_t, boundary_t, MAXORDER> model_t;
typedef Mesh<double, MAXORDER> mesh_t;
typedef LoveMatrices<double, MAXORDER> lovesolver_t;
typedef RayleighMatrices<double, MAXORDER> rayleighsolver_t;


class LeastSquaresIterator {
public:

  LeastSquaresIterator()
  {
  }

  static void copy(const model_t &model, Spec1DMatrix<double> &model_v, Spec1DMatrix<int> &model_mask)
  {
    int k = 0;
    for (auto &c : model.cells) {
      for (int j = 0; j < 4; j ++) {
	for (int i = 0; i <= (int)c.order[j]; i ++) {
	  model_v(k, 0) = c.nodes[i][j];
	  model_mask(k, 0) = j;
	  k ++;
	}
      }
    }
    for (int i = 0; i < 4; i ++) {
      model_v(k, 0) = model.boundary.parameters[i];
      model_mask(k, 0) = i;
      k ++;
    }
  }

  static void initialize_Cm(const model_t &model, const double *damping, Spec1DMatrix<double> &Cm)
  {
    int k = 0;
    for (auto &c : model.cells) {
      for (int j = 0; j < 4; j ++) {
	for (int i = 0; i <= (int)c.order[j]; i ++) {
	  Cm(k, 0) = damping[j] * damping[j];
	  k ++;
	}
      }
    }
    for (int i = 0; i < 4; i ++) {
      Cm(k, 0) = damping[i] * damping[i];
      k ++;
    }
  }

  static void copy(const Spec1DMatrix<double> &model_v, model_t &model)
  {
    int k = 0;
    for (auto &c : model.cells) {
      for (int j = 0; j < 4; j ++) {
	for (int i = 0; i <= (int)c.order[j]; i ++) {
	  c.nodes[i][j] = model_v(k, 0);
	  k ++;
	}
      }
    }
    for (int i = 0; i < 4; i ++) {
      model.boundary.parameters[i] = model_v(k, 0);
      k ++;
    }
  }

  static bool validate(const Spec1DMatrix<double> &model_v,
		       const Spec1DMatrix<int> &model_mask,
		       double *model_vmin,
		       double *model_vmax)
  {
    for (int i = 0; i < model_v.rows(); i ++) {

      int j = model_mask(i, 0);
      if (model_v(i, 0) < model_vmin[j] ||
	  model_v(i, 0) > model_vmax[j]) {
	return false;
      }
    }

    return true;
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
			   Spec1DMatrix<double> &proposed_model) = 0;

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
				Spec1DMatrix<double> &proposed_model) = 0;
				
};


#endif // common_hpp
