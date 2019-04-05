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

double love_phase_compute(DispersionData &data,
			   int i,
			   model_t &model,
			   mesh_t &mesh,
			   lovesolver_t &love,
			   Spec1DMatrix<double> &dkdp,
			   Spec1DMatrix<double> &dUdp,
			   double threshold,
			   int order,
			   int highorder,
			   int boundaryorder,
			   double scale, 
			   double &weight) // out dGm/dm weight for d likelihood/dm)
{
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
  
  love.recompute(mesh, boundaryorder, scale);
  
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

  data.predicted_k[i] = k;
  data.predicted_phase[i] = c_pred;
  data.predicted_group[i] = U_pred;
  
  weight = -(omega/(k*k));
  
  return k;
}

void spline_fill(int offset, int i0, int i1,
		 DispersionData &data,
		 Spec1DMatrix<double> &G,
		 Spec1DMatrix<double> &Gk,
		 Spec1DMatrix<double> &GU)
{
  double omega0 = data.freq[i0] * 2.0 * M_PI;
  double omega1 = data.freq[i1] * 2.0 * M_PI;
  
  for (int i = i0 + 1; i < i1; i ++) {

    double omega = data.freq[i] * 2.0 * M_PI;
    double affine = (omega1 - omega0);
    double t = (omega - omega0)/affine;
    double t2 = t*t;
    double t3 = t2*t;

    /*
     * Interpolate k/c
     */
    double A = 2.0*t3 - 3.0*t2 + 1.0;
    double B = t3 - 2.0*t2 + t;
    double C = -2.0*t3 + 3.0*t2;
    double D = t3 - t2;

    data.predicted_k[i] = A*data.predicted_k[i0] + affine*B/data.predicted_group[i0] +
      C*data.predicted_k[i1] + affine*D/data.predicted_group[i1];

    data.predicted_phase[i] = omega/data.predicted_k[i];

    /*
     * Interpolate Group
     */
    double dA = 6.0*t2 - 6.0*t;
    double dB = 3.0*t2 - 4.0*t + 1.0;
    double dC = -6.0*t2 + 6.0*t;
    double dD = 3.0*t2 - 2.0*t;
    
    data.predicted_group[i] = 1.0/(dA * data.predicted_k[i0]/affine +
				   dB/data.predicted_group[i0] +
				   dC * data.predicted_k[i1]/affine +
				   dD/data.predicted_group[i1]);

    double u20 = data.predicted_group[i0] * data.predicted_group[i0];
    double u21 = data.predicted_group[i1] * data.predicted_group[i1];

    double weight = -(omega)/(data.predicted_k[i]*data.predicted_k[i]);

    /*
     * Compute Interpolated Jacobians
     */
    int datai = i - offset;
    for (int j = 0; j < G.cols(); j ++) {
      Gk(datai, j) =
	A * Gk(i0 - offset, j) -
	affine * B * GU(i0 - offset, j)/u20 +
	C * Gk(i1 - offset, j) -
	affine * D * GU(i1 - offset, j)/u21;
      
      G(datai, j) = weight * Gk(datai, j);
    }
  }
}


double likelihood_love_spline(DispersionData &data,
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
			      Spec1DMatrix<double> &Gk,
			      Spec1DMatrix<double> &GU,
			      Spec1DMatrix<double> &residual,
			      Spec1DMatrix<double> &Cd,
			      double threshold,
			      int order,
			      int highorder,
			      int boundaryorder,
			      double scale,
			      int skip)
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
    size_t ndata = data.flast - data.ffirst + 1;


    //
    // Compue every skip'th frequency
    //
    int i;
    for (i = data.flast; i >= data.ffirst; i -= skip) {

      double weight;
      double k = love_phase_compute(data,
				    i,
				    model,
				    mesh,
				    love,
				    dkdp,
				    dUdp,
				    threshold,
				    order,
				    highorder,
				    boundaryorder,
				    autoscale, 
				    weight);

      if (k == 0.0) {
	return 0.0;
      }

      if (first) {
	G.resize(ndata, dkdp.rows());
	Gk.resize(ndata, dkdp.rows());
	GU.resize(ndata, dkdp.rows());
	residual.resize(ndata, 1);
	Cd.resize(ndata, 1);
      
        dLdp.resize(dkdp.rows(), 1);
        dLdp.setZero();
        first = false;
      }

      int idata = i - data.ffirst;
      
      // residual(idata, 0) = err;
      // Cd(idata, 0) = denom;
      
      for (int j = 0; j < dkdp.rows(); j ++) {

	G(idata, j) = weight * dkdp(j, 0);

	Gk(idata, j) = dkdp(j, 0);
	GU(idata, j) = dUdp(j, 0);
      }
      
      // like += L;

      //
      // Update autoscale
      //
      double vs2 = mesh.boundary.L/mesh.boundary.rho;
      double omega = data.freq[i] * 2.0 * M_PI;
      double disc = k*k - omega*omega/vs2;
      if (disc > 0.0) {
        autoscale = sqrt(disc);
      }
    }

    //
    // Compute first frequency if it was missed in the above loop
    //
    if (i < data.ffirst) {
      i = data.ffirst;
      
      double weight;
      double k = love_phase_compute(data,
				    i,
				    model,
				    mesh,
				    love,
				    dkdp,
				    dUdp,
				    threshold,
				    order,
				    highorder,
				    boundaryorder,
				    autoscale, 
				    weight);

      if (k == 0.0) {
	return 0.0;
      }
	
      for (int j = 0; j < dkdp.rows(); j ++) {

	G(0, j) = weight * dkdp(j, 0);

	Gk(0, j) = dkdp(j, 0);
	GU(0, j) = dUdp(j, 0);
      }
    }

    
    //
    // Fill in missing values of G matrix and predictions with cubic spline
    //
    for (i = data.flast; i >= data.ffirst; i -= skip) {

      int j = i - skip;
      if (j < data.ffirst) {
	j = data.ffirst;
      }

      spline_fill(data.ffirst, j, i, data, G, Gk, GU);
    }

	
    //
    // Compute actual predictions and likelihood dLdp etc from filled G matrix
    //
    for (i = data.flast; i >= data.ffirst; i --) {

      double t_pred = data.predicted_phase[i];
      double err = t_pred - data.target_phase[i];
      double denom = data.target_error[i] * data.target_error[i];
      double L = err*err/(2.0 * denom);
      double normed_residual = err/denom;

      int datai = i - data.ffirst;
      
      if (i == data.flast) {

	for (int j = 0; j < dkdp.rows(); j ++) {
	  dLdp(j, 0) = normed_residual * G(datai, j);
	}
	
      } else {

	for (int j = 0; j < dkdp.rows(); j ++) {
	  dLdp(j, 0) += normed_residual * G(datai, j);
	}
	
      }

      like += L;
    }
  }

  return like;
    
}


double likelihood_love(DispersionData &data,
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
      
      double err = c_pred - data.target_phase[i];
      // if (first) {
      // printf("%15.9f %15.9f %15.9f %15.9f\n", data.freq[i], t_pred, data.time_max[i], err);
      // }
      
      double denom = data.target_error[i] * data.target_error[i];
      
      double L = err*err/(2.0 * denom);

      double weight = -(omega/(k*k));
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

double rayleigh_phase_compute(DispersionData &data,
			      int i,
			      model_t &model,
			      mesh_t &mesh,
			      rayleighsolver_t &rayleigh,
			      Spec1DMatrix<double> &dkdp,
			      Spec1DMatrix<double> &dUdp,
			      double threshold,
			      int order,
			      int highorder,
			      int boundaryorder,
			      double scale, 
			      double &weight) // out dGm/dm weight for d likelihood/dm)
{
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
  
  rayleigh.recompute(mesh, boundaryorder, scale, scale);
  
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
  
  data.predicted_k[i] = k;
  data.predicted_phase[i] = c_pred;
  data.predicted_group[i] = U_pred;
  
  weight = -(omega/(k*k));
  
  return k;
}

double likelihood_rayleigh_spline(DispersionData &data,
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
				  Spec1DMatrix<double> &Gk,
				  Spec1DMatrix<double> &GU,
				  Spec1DMatrix<double> &residual,
				  Spec1DMatrix<double> &Cd,
				  double threshold,
				  int order,
				  int highorder,
				  int boundaryorder,
				  double scale,
				  int skip)
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
    size_t ndata = data.flast - data.ffirst + 1;


    //
    // Compue every skip'th frequency
    //
    int i;
    for (i = data.flast; i >= data.ffirst; i -= skip) {

      double weight;
      double k = rayleigh_phase_compute(data,
					i,
					model,
					mesh,
					rayleigh,
					dkdp,
					dUdp,
					threshold,
					order,
					highorder,
					boundaryorder,
					autoscale, 
					weight);

      if (k == 0.0) {
	return 0.0;
      }

      if (first) {
	G.resize(ndata, dkdp.rows());
	Gk.resize(ndata, dkdp.rows());
	GU.resize(ndata, dkdp.rows());
	residual.resize(ndata, 1);
	Cd.resize(ndata, 1);
      
        dLdp.resize(dkdp.rows(), 1);
        dLdp.setZero();
        first = false;
      }

      int idata = i - data.ffirst;
      
      // residual(idata, 0) = err;
      // Cd(idata, 0) = denom;
      
      for (int j = 0; j < dkdp.rows(); j ++) {

	G(idata, j) = weight * dkdp(j, 0);

	Gk(idata, j) = dkdp(j, 0);
	GU(idata, j) = dUdp(j, 0);
      }
      
      // like += L;

      //
      // Update autoscale
      //
      //
      // Update autoscale: use the vp scale 
      //
      double omega = data.freq[i] * 2.0 * M_PI;
      double vp2 = mesh.boundary.A/mesh.boundary.rho;
      double disc = k*k - omega*omega/vp2;
      if (disc > 0.0) {
        autoscale = sqrt(disc);
      }
    }

    //
    // Compute first frequency if it was missed in the above loop
    //
    if (i < data.ffirst) {
      i = data.ffirst;
      
      double weight;
      double k = rayleigh_phase_compute(data,
					i,
					model,
					mesh,
					rayleigh,
					dkdp,
					dUdp,
					threshold,
					order,
					highorder,
					boundaryorder,
					autoscale, 
					weight);
      
      if (k == 0.0) {
	return 0.0;
      }
	
      for (int j = 0; j < dkdp.rows(); j ++) {

	G(0, j) = weight * dkdp(j, 0);

	Gk(0, j) = dkdp(j, 0);
	GU(0, j) = dUdp(j, 0);
      }
    }

    
    //
    // Fill in missing values of G matrix and predictions with cubic spline
    //
    for (i = data.flast; i >= data.ffirst; i -= skip) {

      int j = i - skip;
      if (j < data.ffirst) {
	j = data.ffirst;
      }

      spline_fill(data.ffirst, j, i, data, G, Gk, GU);
    }

	
    //
    // Compute actual predictions and likelihood dLdp etc from filled G matrix
    //
    for (i = data.flast; i >= data.ffirst; i --) {

      double t_pred = data.predicted_phase[i];
      double err = t_pred - data.target_phase[i];
      double denom = data.target_error[i] * data.target_error[i];
      double L = err*err/(2.0 * denom);
      double normed_residual = err/denom;

      int datai = i - data.ffirst;
      
      if (i == data.flast) {

	for (int j = 0; j < dkdp.rows(); j ++) {
	  dLdp(j, 0) = normed_residual * G(datai, j);
	}
	
      } else {

	for (int j = 0; j < dkdp.rows(); j ++) {
	  dLdp(j, 0) += normed_residual * G(datai, j);
	}
	
      }

      like += L;
    }
  }

  return like;
    
}

double likelihood_rayleigh(DispersionData &data,
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
	for (int j = 0; j < dUdp.rows(); j ++) {
	  dUdp(j, 0) = -dUdp(j, 0);
	}
      }

      data.predicted_phase[i] = c_pred;
      data.predicted_group[i] = U_pred;
      
      double err = c_pred - data.target_phase[i];
      // if (first) {
      // printf("%15.9f %15.9f %15.9f %15.9f\n", data.freq[i], t_pred, data.time_max[i], err);
      // }
      
      double denom = data.target_error[i] * data.target_error[i];

      double L = err*err/(2.0 * denom);
      
      double weight = -omega/(k*k);
      double normed_residual = err/denom;
      
      if (first) {
        dLdp.resize(dkdp.rows(), 1);
        dLdp.setZero();

	G.resize(ndata, dkdp.rows());
	G.setZero();

	Cd.resize(ndata, 1);

	residual.resize(ndata, 1);
	
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


