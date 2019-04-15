//
//    AkiEstimate : A method for the joint estimation of Love and Rayleigh surface wave
//    dispersion from ambient noise cross-correlations.
//
//      Hawkins R. and Sambridge M., "An adjoint technique for estimation of interstation phase
//    and group dispersion from ambient noise cross-correlations", BSSA, 2019
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
#ifndef simple_hpp
#define simple_hpp

class SimpleStep : public LeastSquaresIterator {
public:

  static constexpr double REFERENCE_LEVELS[4] = {
    2500.0, 
    3000.0,
    1.0,
    1.7
  };
  

  
  SimpleStep()
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
    //
    // First compute max gradient for each component
    //
    double maxgradient[4] = {0.0, 0.0, 0.0, 0.0};

    int rows = current_model.rows();

    for (int i = 0; i < rows; i ++) {
      double g = fabs(dLdp(i, 0));
      int j = model_mask(i, 0);
      if (g > maxgradient[j]) {
	maxgradient[j] = g;
      }
    }

    printf("Max gradient  : %16.9e %16.9e %16.9e %16.9e\n", maxgradient[0], maxgradient[1], maxgradient[2], maxgradient[3]);

    double gradient_scale[4];
    for (int i = 0; i < 4; i ++) {
      if (maxgradient[i] > 0.0) {
	gradient_scale[i] = epsilon/100.0 * REFERENCE_LEVELS[i]/maxgradient[i];
      } else {
	gradient_scale[i] = 0.0;
      }
    }
    printf("Gradient Scale: %16.9e %16.9e %16.9e %16.9e\n",
	   gradient_scale[0], gradient_scale[1], gradient_scale[2], gradient_scale[3]);
    
    for (int i = 0; i < rows; i ++) {
      int j = model_mask(i, 0);
      double dp = dLdp(i, 0) * gradient_scale[j];

      proposed_model(i, 0) = current_model(i, 0) - dp;
      if (isnan(proposed_model(i, 0))) {
	fprintf(stderr, "error: new model parameter is NAN %16.9e %16.9e\n", dp, current_model(i, 0));
	return false;
      }
    }

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
    //
    // First compute max gradient for each component
    //
    double maxgradient[4] = {0.0, 0.0, 0.0, 0.0};

    int rows = current_model.rows();

    for (int i = 0; i < rows; i ++) {
      double g = fabs(dLdp(i, 0));
      int j = model_mask(i, 0);
      if (g > maxgradient[j]) {
	maxgradient[j] = g;
      }
    }

    printf("Max gradient  : %16.9e %16.9e %16.9e %16.9e\n", maxgradient[0], maxgradient[1], maxgradient[2], maxgradient[3]);

    double gradient_scale[4];
    for (int i = 0; i < 4; i ++) {
      if (maxgradient[i] > 0.0) {
	gradient_scale[i] = epsilon/100.0 * REFERENCE_LEVELS[i]/maxgradient[i];
      } else {
	gradient_scale[i] = 0.0;
      }
    }
    printf("Gradient Scale: %16.9e %16.9e %16.9e %16.9e\n",
	   gradient_scale[0], gradient_scale[1], gradient_scale[2], gradient_scale[3]);
    
    for (int i = 0; i < rows; i ++) {
      int j = model_mask(i, 0);
      double dp = dLdp(i, 0) * gradient_scale[j];

      proposed_model(i, 0) = current_model(i, 0) - dp;
    }

    return true;
  }
  
};

#endif // simple_hpp
