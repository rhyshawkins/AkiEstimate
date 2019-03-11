#ifndef likelihood_hpp
#define likelihood_hpp

#include "common.hpp"
#include "dispersion.hpp"
#include "reference.hpp"

void likelihood_damping(model_t &model,
			model_t &reference,
			const double *damping,
			Spec1DMatrix<double> &dLdp,
			double &like)
{
  //
  // Cells penalties
  //
  int offset = 0;
  int ci = 0;
  for (auto &c : reference.cells) {
    
    for (int j = 0; j < 4; j ++) {
      
      if (damping[j] > 0.0) {
	
	double denom = damping[j] * damping[j];
	
	for (int i = 0; i <= (int)c.order[j]; i ++) {
	  
	  double dm = (model.cells[ci].nodes[i][j] - c.nodes[i][j]);

	  // printf("%d: %10.6f %10.6f : %10.6f/%10.6f\n", j,
	  // 	 model.cells[ci].nodes[i][j], c.nodes[i][j],
	  // 	 dm, damping[j]);
	  
	  like += (dm*dm)/(2.0 * denom);
	  
	  dLdp(offset + i, 0) -= dm/denom;

	}
	
      }

      offset += c.order[j] + 1;

    }
    
    ci ++;
  }

  //
  // Halfspace
  //

  for (int i = 0; i < 4; i ++) {

    if (damping[i] > 0.0) {
      double dm =  model.boundary.parameters[i] - reference.boundary.parameters[i];
      double denom = damping[i] * damping[i];
      like += (dm*dm)/(2.0 * denom);
      dLdp(offset + i, 0) -= dm/denom;
    }
  }
}

bool love_jacobian(DispersionData &data,
		   model_t &model,
		   model_t &reference,
		   const double *damping,
		   mesh_t &mesh,
		   lovesolver_t &love,
		   Spec1DMatrix<double> &dkdp,
		   Spec1DMatrix<double> &dUdp,
		   Spec1DMatrix<double> &Jc,
		   Spec1DMatrix<double> &JU,
		   double threshold,
		   int order,
		   int highorder,
		   int boundaryorder,
		   double scale,
		   double frequency_thin)
{
  double autoscale = scale;
  bool first = true;

  double last_freq = -1.0;
  
  //
  // Count actual no. data
  //
  size_t ndata = 0;
  for (int i = data.flast; i >= data.ffirst; i --) {
    
    if (frequency_thin > 0.0 && last_freq > 0.0) {
      if (last_freq - data.freq[i] < frequency_thin) {
	// printf("%15.9f skipped\n", data.freq[i]);
	continue;
      }
    }
    
    last_freq = data.freq[i];
    ndata ++;
  }
  
  last_freq = -1.0;
  int data_i = 0;
  
  for (int i = data.flast; i >= data.ffirst; i --) {
    
    if (frequency_thin > 0.0 && last_freq > 0.0) {
      if (last_freq - data.freq[i] < frequency_thin) {
	  // printf("%15.9f skipped\n", data.freq[i]);
	continue;
      }
    }
    
    last_freq = data.freq[i];
    
    if (threshold <= 0.0) {
      
      model.project_gradient(mesh, order);
      
    } else {
      printf("Unimplemented\n");
      model.project_threshold(mesh,
			      threshold,
			      highorder,
			      0.0,
			      order);
      
    }
    
    love.recompute(mesh, boundaryorder, autoscale);
    
    double normA, normB, normC;
    double omega = data.freq[i] * 2.0 * M_PI;
    double k = love.solve_fundamental_gradient_sep(mesh,
						   boundaryorder,
						   omega,
						   dkdp,
						   dUdp,
						   normA,
						   normB,
						   normC);
    if (k <= 0.0) {
      fprintf(stderr, "error: failed to compute wave number (%f, %f %f %f)\n",
	      omega,
	      love.A(0, 0),
	      love.B(0, 0),
	      love.C(0, 0));
      exit(-1);
    }
    
    if (first) {
      Jc.resize(ndata, dkdp.rows());
      JU.resize(ndata, dkdp.rows());
      first = false;
    }
    
    for (int j = 0; j < dkdp.rows(); j ++) {
      Jc(data_i, j) = -omega/(k*k) * dkdp(j, 0);
      JU(data_i, j) = dUdp(j, 0);
    }
      
    data_i ++;

    //
    // Update autoscale
    //
    double vs2 = mesh.boundary.L/mesh.boundary.rho;
    double disc = k*k - omega*omega/vs2;
    if (disc > 0.0) {
      autoscale = sqrt(disc);
    }
  }

  return true;
}


double likelihood_love_bessel(DispersionData &data,
			      model_t &model,
			      model_t &reference,
			      const double *damping,
			      bool posterior,
			      mesh_t &mesh,
			      lovesolver_t &love,
			      Spec1DMatrix<double> &dkdp,
			      Spec1DMatrix<double> &dUdp,
			      Spec1DMatrix<double> &dLdp,
			      Spec1DMatrix<double> &G,
			      Spec1DMatrix<double> &residual,
			      Spec1DMatrix<double> &Cd,
			      double threshold,
			      int order,
			      int highorder,
			      int boundaryorder,
			      double scale,
			      double frequency_thin)
{
  double autoscale = scale;
  bool first = true;
  double like = 0.0;

  if (posterior) {

    if (dLdp.rows() == 0) {

      int count = 1;
      for (auto &c : model.cells) {
        count +=
	  (c.order[0] + 1) +
	  (c.order[1] + 1) +
	  (c.order[2] + 1) +
	  (c.order[3] + 1);
      }

      count += 4; // For halfspace
      
      dLdp.resize(count, 1);
    }

    dLdp.setZero();
    
  } else {

    double last_freq = -1.0;
    
    //
    // Count actual no. data
    //
    size_t ndata = 0;
    for (int i = data.flast; i >= data.ffirst; i --) {

      if (frequency_thin > 0.0 && last_freq > 0.0) {
	if (last_freq - data.freq[i] < frequency_thin) {
	  // printf("%15.9f skipped\n", data.freq[i]);
	  continue;
	}
      }

      last_freq = data.freq[i];
      ndata ++;
    }

    last_freq = -1.0;
    int data_i = 0;

    for (int i = data.flast; i >= data.ffirst; i --) {

      if (frequency_thin > 0.0 && last_freq > 0.0) {
	if (last_freq - data.freq[i] < frequency_thin) {
	  // printf("%15.9f skipped\n", data.freq[i]);
	  continue;
	}
      }

      last_freq = data.freq[i];
      
      if (threshold <= 0.0) {
	
        model.project_gradient(mesh, order);
        
      } else {
        printf("Unimplemented\n");
        model.project_threshold(mesh,
                                threshold,
                                highorder,
                                0.0,
                                order);
	
      }
      
      love.recompute(mesh, boundaryorder, autoscale);
      
      double normA, normB, normC;
      double omega = data.freq[i] * 2.0 * M_PI;
      double k = love.solve_fundamental_gradient_sep(mesh,
						     boundaryorder,
						     omega,
						     dkdp,
						     dUdp,
						     normA,
						     normB,
						     normC);
      if (k <= 0.0) {
        fprintf(stderr, "error: failed to compute wave number (%f, %f %f %f)\n",
                omega,
                love.A(0, 0),
                love.B(0, 0),
                love.C(0, 0));
        exit(-1);
      }
      
      double c_pred = omega/k;
      double U_pred = (k * normB)/(omega * normA);

      data.predicted_phase[i] = c_pred;
      data.predicted_group[i] = U_pred;

      double pJ0 = gsl_sf_bessel_J0(k * data.distkm * 1.0e3);
      double pJ1 = gsl_sf_bessel_J1(k * data.distkm * 1.0e3);
      double pY0 = gsl_sf_bessel_Y0(k * data.distkm * 1.0e3);
      double pY1 = gsl_sf_bessel_Y1(k * data.distkm * 1.0e3);
      
      data.predicted_bessel[i] = pJ0;

      double benv2 = pJ0*pJ0 + pY0*pY0;
      double benv = sqrt(benv2);
      double s = data.predicted_envelope[i]/benv;
      data.predicted_realspec[i] = s * data.predicted_bessel[i];
      
      double t_pred = data.predicted_realspec[i];

      double err = t_pred - data.ncfreal[i];
      // if (first) {
      // printf("%15.9f %15.9f %15.9f %15.9f\n", data.freq[i], t_pred, data.time_max[i], err);
      // }
      
      double denom = data.noise_sigma * data.noise_sigma;
      double L = err*err/(2.0 * denom);

      double dpbdk = data.predicted_envelope[i] *
	(- (pJ1 * data.distkm * 1.0e3)/benv
	 - (pJ0*(pJ0*pJ1 + pY0*pY1) * data.distkm * 1.0e3)/(benv2*benv));
      
      double weight = dpbdk;
      double normed_residual = err/denom;
      
      if (first) {
	G.resize(ndata, dkdp.rows());
	residual.resize(ndata, 1);
	Cd.resize(ndata, 1);
      
        dLdp.resize(dkdp.rows(), 1);
        dLdp.setZero();
        first = false;
      }

      residual(data_i, 0) = err;
      Cd(data_i, 0) = denom;
      
      for (int j = 0; j < dkdp.rows(); j ++) {
        dLdp(j, 0) += weight * normed_residual * dkdp(j, 0);

	G(data_i, j) = weight * dkdp(j, 0);
      }
      
      like += L;
      data_i ++;

      //
      // Update autoscale
      //
      double vs2 = mesh.boundary.L/mesh.boundary.rho;
      double disc = k*k - omega*omega/vs2;
      if (disc > 0.0) {
        autoscale = sqrt(disc);
      }
    }
  }

  likelihood_damping(model,
		     reference,
		     damping,
		     dLdp,
		     like);

  return like;
}

bool rayleigh_jacobian(DispersionData &data,
		       model_t &model,
		       model_t &reference,
		       const double *damping,
		       mesh_t &mesh,
		       rayleighsolver_t &rayleigh,
		       Spec1DMatrix<double> &dkdp,
		       Spec1DMatrix<double> &dUdp,
		       Spec1DMatrix<double> &Jc,
		       Spec1DMatrix<double> &JU,
		       double threshold,
		       int order,
		       int highorder,
		       int boundaryorder,
		       double scale,
		       double frequency_thin)
{
  double autoscale = scale;
  bool first = true;

  double last_freq = -1.0;
  
  //
  // Count actual no. data
  //
  size_t ndata = 0;
  for (int i = data.flast; i >= data.ffirst; i --) {
    
    if (frequency_thin > 0.0 && last_freq > 0.0) {
      if (last_freq - data.freq[i] < frequency_thin) {
	// printf("%15.9f skipped\n", data.freq[i]);
	continue;
      }
    }
    
    last_freq = data.freq[i];
    ndata ++;
  }
  
  last_freq = -1.0;
  int data_i = 0;
  
  for (int i = data.flast; i >= data.ffirst; i --) {
    
    if (frequency_thin > 0.0 && last_freq > 0.0) {
      if (last_freq - data.freq[i] < frequency_thin) {
	// printf("%15.9f skipped\n", data.freq[i]);
	continue;
      }
    }
    
    last_freq = data.freq[i];
    
    if (threshold <= 0.0) {
      
      model.project_gradient(mesh, order);
      
    } else {
      printf("Unimplemented\n");
      model.project_threshold(mesh,
			      threshold,
			      highorder,
			      0.0,
			      order);
      
    }
    
    rayleigh.recompute(mesh, boundaryorder, autoscale, autoscale);
    
    double normA, normB, normC, normD;
    double omega = data.freq[i] * 2.0 * M_PI;
    Spec1DMatrix<double> dgxdv;
    Spec1DMatrix<double> dgzdv;
    Spec1DMatrix<double> dGvdp;
    double _eH;
    double _eV;
    
    dgxdv.resize(rayleigh.size, 1);
    dgxdv.setZero();
    dgxdv(0, 0) = 1.0;
    
    dgzdv.resize(rayleigh.size, 1);
    dgzdv.setZero();
    dgzdv(0, 0) = 1.0;
    
    double k = rayleigh.solve_fundamental_gradient_generic(mesh,
							   boundaryorder,
							   omega,
							   dgxdv,
							   dgzdv,
							   dkdp,
							   dUdp,
							   normA,
							   normB,
							   normC,
							   normD,
							   _eH,
							   _eV,
							   dGvdp);
    if (k == 0.0) {
      fprintf(stderr, "error: failed to compute wave number (%f, %f %f %f)\n",
	      omega,
	      rayleigh.Ax(0, 0),
	      rayleigh.Bx(0, 0),
	      rayleigh.Cx(0, 0));
      exit(-1);
    }
    
    double c_pred = omega/k;
    double U_pred = (2.0*normB*k + normC)/(2.0*omega*normA);
    if (k < 0.0) {
      U_pred = -U_pred;
      c_pred = -c_pred;
      for (int j = 0; j < dUdp.rows(); j ++) {
	dUdp(j, 0) = -dUdp(j, 0);
      }
      for (int j = 0; j < dkdp.rows(); j ++) {
	dkdp(j, 0) = -dkdp(j, 0);
      }
      
      k = -k;
    }
    
    if (first) {
      Jc.resize(ndata, dkdp.rows());
      JU.resize(ndata, dkdp.rows());
      first = false;
    }
      
    for (int j = 0; j < dkdp.rows(); j ++) {
      Jc(data_i, j) = -omega/(k*k) * dkdp(j, 0);
      JU(data_i, j) = dUdp(j, 0);
    }
    
    data_i ++;
      
    //
    // Update autoscale: use the vp scale 
    //
    double vp2 = mesh.boundary.A/mesh.boundary.rho;
    double disc = k*k - omega*omega/vp2;
    if (disc > 0.0) {
      autoscale = sqrt(disc);
    }
  }

  return true;
}

double likelihood_rayleigh_bessel(DispersionData &data,
				  model_t &model,
				  model_t &reference,
				  const double *damping,
				  bool posterior,
				  mesh_t &mesh,
				  rayleighsolver_t &rayleigh,
				  Spec1DMatrix<double> &dkdp,
				  Spec1DMatrix<double> &dUdp,
				  Spec1DMatrix<double> &dLdp,
				  Spec1DMatrix<double> &G,
				  Spec1DMatrix<double> &residual,
				  Spec1DMatrix<double> &Cd,
				  double threshold,
				  int order,
				  int highorder,
				  int boundaryorder,
				  double scale,
				  double frequency_thin)
{
  double autoscale = scale;
  bool first = true;
  double like = 0.0;

  if (posterior) {

    if (dLdp.rows() == 0) {

      int count = 1;
      for (auto &c : model.cells) {
        count +=
	  (c.order[0] + 1) +
	  (c.order[1] + 1) +
	  (c.order[2] + 1) +
	  (c.order[3] + 1);
      }

      count += 4; // For halfspace
      
      dLdp.resize(count, 1);
    }

    dLdp.setZero();
    
  } else {

    double last_freq = -1.0;
    
    //
    // Count actual no. data
    //
    size_t ndata = 0;
    for (int i = data.flast; i >= data.ffirst; i --) {

      if (frequency_thin > 0.0 && last_freq > 0.0) {
	if (last_freq - data.freq[i] < frequency_thin) {
	  // printf("%15.9f skipped\n", data.freq[i]);
	  continue;
	}
      }

      last_freq = data.freq[i];
      ndata ++;
    }

    last_freq = -1.0;
    int data_i = 0;

    for (int i = data.flast; i >= data.ffirst; i --) {

      if (frequency_thin > 0.0 && last_freq > 0.0) {
	if (last_freq - data.freq[i] < frequency_thin) {
	  // printf("%15.9f skipped\n", data.freq[i]);
	  continue;
	}
      }

      last_freq = data.freq[i];
      
      if (threshold <= 0.0) {
	
        model.project_gradient(mesh, order);
        
      } else {
        printf("Unimplemented\n");
        model.project_threshold(mesh,
                                threshold,
                                highorder,
                                0.0,
                                order);
	
      }
      
      rayleigh.recompute(mesh, boundaryorder, autoscale, autoscale);
      
      double normA, normB, normC, normD;
      double omega = data.freq[i] * 2.0 * M_PI;
      Spec1DMatrix<double> dgxdv;
      Spec1DMatrix<double> dgzdv;
      Spec1DMatrix<double> dGvdp;
      double _eH;
      double _eV;

      dgxdv.resize(rayleigh.size, 1);
      dgxdv.setZero();
      dgxdv(0, 0) = 1.0;
      
      dgzdv.resize(rayleigh.size, 1);
      dgzdv.setZero();
      dgzdv(0, 0) = 1.0;
      
      double k = rayleigh.solve_fundamental_gradient_generic(mesh,
							     boundaryorder,
							     omega,
							     dgxdv,
							     dgzdv,
							     dkdp,
							     dUdp,
							     normA,
							     normB,
							     normC,
							     normD,
							     _eH,
							     _eV,
							     dGvdp);
      if (k == 0.0) {
        fprintf(stderr, "error: failed to compute wave number (%f, %f %f %f)\n",
                omega,
                rayleigh.Ax(0, 0),
                rayleigh.Bx(0, 0),
                rayleigh.Cx(0, 0));
        exit(-1);
      }
      
      double c_pred = omega/k;
      double U_pred = (2.0*normB*k + normC)/(2.0*omega*normA);
      if (k < 0.0) {
	U_pred = -U_pred;
	c_pred = -c_pred;
	for (int j = 0; j < dUdp.rows(); j ++) {
	  dUdp(j, 0) = -dUdp(j, 0);
	}
	for (int j = 0; j < dkdp.rows(); j ++) {
	  dkdp(j, 0) = -dkdp(j, 0);
	}

	k = -k;
      }
      
      data.predicted_phase[i] = c_pred;
      data.predicted_group[i] = U_pred;

      double pJ0 = gsl_sf_bessel_J0(k * data.distkm * 1.0e3);
      double pJ1 = gsl_sf_bessel_J1(k * data.distkm * 1.0e3);
      double pY0 = gsl_sf_bessel_Y0(k * data.distkm * 1.0e3);
      double pY1 = gsl_sf_bessel_Y1(k * data.distkm * 1.0e3);
      
      data.predicted_bessel[i] = pJ0;

      double benv2 = pJ0*pJ0 + pY0*pY0;
      double benv = sqrt(benv2);
      double s = data.predicted_envelope[i]/benv;
      data.predicted_realspec[i] = s * data.predicted_bessel[i];
      
      double t_pred = data.predicted_realspec[i];

      double err = t_pred - data.ncfreal[i];
      
      double denom = data.noise_sigma * data.noise_sigma;
      double L = err*err/(2.0 * denom);
      
      double dpbdk = data.predicted_envelope[i] *
	(- (pJ1 * data.distkm * 1.0e3)/benv
	 - (pJ0*(pJ0*pJ1 + pY0*pY1) * data.distkm * 1.0e3)/(benv2*benv));
      
      double weight = dpbdk;
      double normed_residual = err/denom;
      
      if (first) {
	G.resize(ndata, dkdp.rows());
	Cd.resize(ndata, 1);
	residual.resize(ndata, 1);

        dLdp.resize(dkdp.rows(), 1);
        dLdp.setZero();
        first = false;
      }
      
      residual(data_i, 0) = err;
      Cd(data_i, 0) = denom;

      for (int j = 0; j < dkdp.rows(); j ++) {
        dLdp(j, 0) += weight * normed_residual * dkdp(j, 0);
	G(data_i, j) = weight * dkdp(j, 0);
      }
      
      like += L;
      data_i ++;
      
      //
      // Update autoscale: use the vp scale 
      //
      double vp2 = mesh.boundary.A/mesh.boundary.rho;
      double disc = k*k - omega*omega/vp2;
      if (disc > 0.0) {
        autoscale = sqrt(disc);
      }
    }
  }

  likelihood_damping(model,
		     reference,
		     damping,
		     dLdp,
		     like);

  return like;
}


#endif // likelihood_hpp


