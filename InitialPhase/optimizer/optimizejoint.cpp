#include <exception>

#include <stdio.h>
#include <getopt.h>

#include "likelihood.hpp"
#include "autosigma.hpp"

#include "simple.hpp"
#include "quasinewton.hpp"

static char short_options[] = "i:c:I:C:r:f:F:R:V:X:S:o:s:p:b:t:P:N:e:QA:WM:T:h";
static struct option long_options[] = {
  {"input-love", required_argument, 0, 'i'},
  {"phase-love", required_argument, 0, 'c'},
  {"input-rayleigh", required_argument, 0, 'I'},
  {"phase-rayleigh", required_argument, 0, 'C'},
  
  {"reference", required_argument, 0, 'r'},

  {"fmin", required_argument, 0, 'f'},
  {"fmax", required_argument, 0, 'F'},

  {"sigma-rho", required_argument, 0, 'R'},
  {"sigma-vs", required_argument, 0, 'V'},
  {"sigma-xi", required_argument, 0, 'X'},
  {"sigma-vpvs", required_argument, 0, 'S'},

  {"output", required_argument, 0, 'o'},

  {"scale", required_argument, 0, 's'},
  {"order", required_argument, 0, 'p'},
  {"boundaryorder", required_argument, 0, 'b'},
  {"threshold", required_argument, 0, 't'},
  {"high-order", required_argument, 0, 'P'},

  {"nsteps", required_argument, 0, 'N'},
  {"epsilon", required_argument, 0, 'e'},

  {"posterior", no_argument, 0, 'Q'},
  
  {"amplitude-threshold", required_argument, 0, 'A'},
  {"write-ftan", no_argument, 0, 'W'},
  {"mode", required_argument, 0, 'M'},

  {"thin", required_argument, 0, 'T'},

  {"help", no_argument, 0, 'h'},
  
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static bool invert(DispersionData &data_love,
		   DispersionData &data_rayleigh,
                   model_t &model,
                   model_t &reference,
                   const double *damping,
                   bool posterior,
                   mesh_t &mesh,
                   lovesolver_t &love,
		   rayleighsolver_t &rayleigh,
                   double threshold,
                   double order,
                   double highorder,
                   double boundaryorder,
                   double scale,
                   double epsilon,
                   int maxiterations,
                   double likelihood_threshold,
		   int mode,
		   int skip);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  char *input_love;
  char *phase_love;
  char *input_rayleigh;
  char *phase_rayleigh;
  
  char *reference_file;
  char *output_file;

  double fmin;
  double fmax;

  double threshold;
  int order;
  int highorder;
  int boundaryorder;
  double scale;

  int maxiterations;
  double epsilon;
  double damping[4];

  bool nodata;

  double sigma;
  double vmin;
  double vmax;
  double max_deltav;

  double amplitude_threshold;

  bool write_ftan;

  int mode;

  int skip;

  //
  // Defaults
  //
  input_love = nullptr;
  phase_love = nullptr;
  
  input_rayleigh = nullptr;
  phase_rayleigh = nullptr;
  
  reference_file = nullptr;
  output_file = nullptr;
  threshold = 0.0;

  order = 5;
  highorder = 5;
  boundaryorder = 5;

  scale = 1.0e-4;
  
  fmin = 1.0/40.0;
  fmax = 1.0/2.0;
  
  maxiterations = 5;
  epsilon = 1.0;

  damping[0] = 0.0;
  damping[1] = 0.0;
  damping[2] = 0.0;
  damping[3] = 0.0;

  nodata = false;

  vmin = 1.0;
  vmax = 5.0;

  max_deltav = 1.0e-3;

  amplitude_threshold = 0.0;

  write_ftan = false;

  mode = 0;

  skip = 0;
  
  //
  // Command line parameters
  //
  option_index = 0;
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input_love = optarg;
      break;

    case 'c':
      phase_love = optarg;
      break;

    case 'I':
      input_rayleigh = optarg;
      break;

    case 'C':
      phase_rayleigh = optarg;
      break;
      
    case 'r':
      reference_file = optarg;
      break;
      
    case 'o':
      output_file = optarg;
      break;

    case 'f':
      fmin = atof(optarg);
      if (fmin <= 0.0) {
	fprintf(stderr, "error: fmin must be greater than 0\n");
	return -1;
      }
      break;

    case 'F':
      fmax = atof(optarg);
      if (fmax <= 0.0) {
	fprintf(stderr, "error: fmax must be greater than 0\n");
	return -1;
      }
      break;

    case 'R':
      damping[0] = atof(optarg);
      if (damping[0] < 0.0) {
	fprintf(stderr, "error: rho std-dev must be 0 or greater\n");
	return -1;
      }
      break;

    case 'V':
      damping[1] = atof(optarg);
      if (damping[1] < 0.0) {
	fprintf(stderr, "error: vs std-dev must be 0 or greater\n");
	return -1;
      }
      break;

    case 'X':
      damping[2] = atof(optarg);
      if (damping[2] < 0.0) {
	fprintf(stderr, "error: xi std-dev must be 0 or greater\n");
	return -1;
      }
      break;

    case 'S':
      damping[3] = atof(optarg);
      if (damping[3] < 0.0) {
	fprintf(stderr, "error: vp/vs std-dev must be 0 or greater\n");
	return -1;
      }
      break;

    case 's':
      scale = atof(optarg);
      if (scale <= 0.0) {
        fprintf(stderr, "error: scale must be positive\n");
        return -1;
      }
      break;

    case 'p':
      order = atoi(optarg);
      if (order < 1) {
        fprintf(stderr, "error: order must be 1 or greater\n");
        return -1;
      }
      break;
      
    case 'b':
      boundaryorder = atoi(optarg);
      if (boundaryorder < 1) {
        fprintf(stderr, "error: boundary order must be 1 or greater\n");
        return -1;
      }
      break;

    case 't':
      threshold = atof(optarg);
      break;

    case 'P':
      highorder = atoi(optarg);
      if (highorder < 1) {
        fprintf(stderr, "error: high order must be 1 or greater\n");
        return -1;
      }
      break;

    case 'N':
      maxiterations = atoi(optarg);
      if (maxiterations < 1) {
        fprintf(stderr, "error: need at least one iteration\n");
        return -1;
      }
      break;

    case 'e':
      epsilon = atof(optarg);
      if (epsilon <= 0.0) {
        fprintf(stderr, "error: epsilon must be positive\n"); 
        return -1;
      }
      break;

    case 'Q':
      nodata = true;
      break;

    case 'A':
      amplitude_threshold = atof(optarg);
      break;

    case 'W':
      write_ftan = true;
      break;
      
    case 'M':
      mode = atoi(optarg);
      if (mode < 0 || mode > 1) {
	fprintf(stderr, "error: mode must be 0 (simple gradient desc.) or 1 (q-newton)\n");
	return -1;
      }
      break;

    case 'T':
      skip = atoi(optarg);
      if (skip < 0) {
	fprintf(stderr, "error: thin must be positive\n");
	return -1;
      }
      break;

    default:
      fprintf(stderr, "unknown option %c\n", c);
    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (input_love == nullptr) {
    fprintf(stderr, "error: missing input love file paramter\n");
    return -1;
  }

  if (phase_love == nullptr) {
    fprintf(stderr, "error: missing input phase love file parameter\n");
    return -1;
  }

  if (input_rayleigh == nullptr) {
    fprintf(stderr, "error: missing input rayleigh file parameter\n");
    return -1;
  }

  if (phase_rayleigh == nullptr) {
    fprintf(stderr, "error: missing input phase rayleigh file parameter\n");
    return -1;
  }
  
  if (reference_file == nullptr) {
    fprintf(stderr, "error: missing reference file parameter\n");
    return -1;
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: missing output file parameter\n");
    return -1;
  }

  DispersionData data_love(fmin, fmax);
  DispersionData data_rayleigh(fmin, fmax);

  if (!data_love.load(input_love)) {
    fprintf(stderr, "error: failed to load love data\n");
    return -1;
  }

  if (!data_love.load_phase(phase_love)) {
    return -1;
  }

  printf("Love Desired range: %10.6f %10.6f\n",
	 data_love.freq[data_love.ffirst],
	 data_love.freq[data_love.flast]);
  data_love.initialise_target();
  printf("Love Actual  range: %10.6f %10.6f\n",
	 data_love.freq[data_love.ffirst],
	 data_love.freq[data_love.flast]);

  if (!data_rayleigh.load(input_rayleigh)) {
    fprintf(stderr, "error: failed to load rayleigh data\n");
    return -1;
  }

  if (!data_rayleigh.load_phase(phase_rayleigh)) {
    return -1;
  }

  printf("Rayleigh Desired range: %10.6f %10.6f\n",
	 data_rayleigh.freq[data_rayleigh.ffirst],
	 data_rayleigh.freq[data_rayleigh.flast]);
  data_rayleigh.initialise_target();
  printf("Rayleigh Actual  range: %10.6f %10.6f\n",
	 data_rayleigh.freq[data_rayleigh.ffirst],
	 data_rayleigh.freq[data_rayleigh.flast]);

  Mesh<double, MAXORDER> mesh;
  MeshAmplitude<double, MAXORDER, MAXORDER> amplitude;

  //
  // Load reference model
  //
  ReferenceModel reference;
  bool promote = true;
  size_t promote_order = order;

  if (!reference.load(reference_file, promote, promote_order)) {
    fprintf(stderr, "error: failed to load model from %s\n", reference_file);
    return -1;
  }

  LoveMatrices<double, MAXORDER, BOUNDARYORDER> love;
  RayleighMatrices<double, MAXORDER, BOUNDARYORDER> rayleigh;
  
  if (!invert(data_love,
	      data_rayleigh,
	      reference.model,
	      reference.reference,
	      damping,
	      nodata,
	      mesh,
	      love,
	      rayleigh,
	      threshold,
	      order,
	      highorder,
	      boundaryorder,
	      scale,
	      epsilon,
	      maxiterations,
	      0.0,
	      mode,
	      skip)) {
    fprintf(stderr, "error: failed to invert\n");
    return -1;
  }

  char filename[1024];
  
  //
  // Save model
  //
  sprintf(filename, "%s.model", output_file);
  if (!reference.model.save(filename)) {
    fprintf(stderr, "error: failed to save model\n");
    return -1;
  }
  
  //
  // Save residuals
  //
  
  
  
  //
  // Save predictions
  //
  sprintf(filename, "%s.pred-love", output_file);
  if (!data_love.save_predictions(filename)) {
    fprintf(stderr, "error: failed to save predictions\n");
    return -1;
  }

  sprintf(filename, "%s.pred-rayleigh", output_file);
  if (!data_rayleigh.save_predictions(filename)) {
    fprintf(stderr, "error: failed to save predictions\n");
    return -1;
  }

  return 0;
}

void usage(const char *pname)
{
  fprintf(stderr,
          "usage: %s [options]\n"
          "where options is one or more of:\n"
          "\n"
          " -i|--input <filename>           Input model (required)\n"
          " -f|--frequency <float>          Frequency\n"
          " -s|--scale <float>              Laguerre scaling (initial)\n"
          "\n"
          " -h|--help                       Show usage information\n"
          "\n",
          pname);
}


static bool invert(DispersionData &data_love,
		   DispersionData &data_rayleigh,
                   model_t &model,
                   model_t &reference,
                   const double *damping,
                   bool posterior,
                   mesh_t &mesh,
                   lovesolver_t &love,
		   rayleighsolver_t &rayleigh,
                   double threshold,
                   double order,
                   double highorder,
                   double boundaryorder,
                   double scale,
                   double _epsilon,
                   int maxiterations,
                   double likelihood_threshold,
		   int mode,
		   int skip)
{
  Spec1DMatrix<double> dkdp_love;
  Spec1DMatrix<double> dUdp_love;
  Spec1DMatrix<double> dLdp_love;
  Spec1DMatrix<double> G_love;
  Spec1DMatrix<double> Gk_love;
  Spec1DMatrix<double> GU_love;
  Spec1DMatrix<double> residuals_love;
  Spec1DMatrix<double> Cd_love;

  Spec1DMatrix<double> old_residuals_love;
  Spec1DMatrix<double> old_residuals_rayleigh;
  Spec1DMatrix<double> old_G_love;
  Spec1DMatrix<double> old_G_rayleigh;
  Spec1DMatrix<double> old_dLdp_love;

  Spec1DMatrix<double> dkdp_rayleigh;
  Spec1DMatrix<double> dUdp_rayleigh;
  Spec1DMatrix<double> dLdp_rayleigh;
  Spec1DMatrix<double> G_rayleigh;
  Spec1DMatrix<double> Gk_rayleigh;
  Spec1DMatrix<double> GU_rayleigh;
  Spec1DMatrix<double> residuals_rayleigh;
  Spec1DMatrix<double> Cd_rayleigh;

  Spec1DMatrix<int> model_mask;
  Spec1DMatrix<double> model_v;
  Spec1DMatrix<double> model_v_proposed;
  Spec1DMatrix<double> model_0;
  Spec1DMatrix<double> Cm;
  
  double frequency_thin = 0.001;

  double epsilon[2];
  LeastSquaresIterator *step[2];

  double PRIOR_MIN[4] = {0.1e3, 0.5e3, 0.5, 1.0};
  double PRIOR_MAX[4] = {8.0e3, 10.0e3, 1.5, 2.5};

  epsilon[0] = _epsilon;
  epsilon[1] = _epsilon/8.0;

  step[0] = new SimpleStep();
  step[1] = new QuasiNewton();
    
  double like_love;
  double like_rayleigh;

  if (skip <= 1) {

    like_love = likelihood_love(data_love,
				model,
				reference,
				damping,
				posterior,
				mesh,
				love,
				dkdp_love,
				dUdp_love,
				dLdp_love,
				G_love,
				residuals_love,
				Cd_love,
				threshold,
				order,
				highorder,
				boundaryorder,
				scale,
				frequency_thin);
    
    like_rayleigh = likelihood_rayleigh(data_rayleigh,
					model,
					reference,
					damping,
					posterior,
					mesh,
					rayleigh,
					dkdp_rayleigh,
					dUdp_rayleigh,
					dLdp_rayleigh,
					G_rayleigh,
					residuals_rayleigh,
					Cd_rayleigh,
					threshold,
					order,
					highorder,
					boundaryorder,
					scale,
					frequency_thin);
  } else {

    like_love = likelihood_love_spline(data_love,
				       model,
				       reference,
				       damping,
				       false,
				       mesh,
				       love,
				       dkdp_love,
				       dUdp_love,
				       dLdp_love,
				       G_love,
				       Gk_love,
				       GU_love,
				       residuals_love,
				       Cd_love,
				       threshold,
				       order,
				       highorder,
				       boundaryorder,
				       scale,
				       skip);
    
    like_rayleigh = likelihood_rayleigh_spline(data_rayleigh,
					       model,
					       reference,
					       damping,
					       false,
					       mesh,
					       rayleigh,
					       dkdp_rayleigh,
					       dUdp_rayleigh,
					       dLdp_rayleigh,
					       G_rayleigh,
					       Gk_rayleigh,
					       GU_rayleigh,
					       residuals_rayleigh,
					       Cd_rayleigh,
					       threshold,
					       order,
					       highorder,
					       boundaryorder,
					       scale,
					       skip);
  }
    
  double like = like_love + like_rayleigh;
  printf("init: %16.9e\n", like);
  double last_like = like;

  //
  // Resize vectors/matrices G will be filled appropriately by likelihood
  // and will be correct size.
  //
  size_t nparam = G_love.cols();

  //
  // Store parameters for perturbing current model
  //
  old_residuals_love = residuals_love;
  old_residuals_rayleigh = residuals_rayleigh;
  old_G_love = G_love;
  old_G_rayleigh = G_rayleigh;

  for (size_t i = 0; i < nparam; i ++) {
    dLdp_love(i, 0) += dLdp_rayleigh(i, 0);
  }

  old_dLdp_love = dLdp_love;

  //
  // Diagonal model covariance matrices
  //
  Cm.resize(nparam, 1);
  LeastSquaresIterator::initialize_Cm(model, damping, Cm);
  
  //
  // Model vectors
  //
  model_mask.resize(nparam, 1);
  model_0.resize(nparam, 1);
  model_v.resize(nparam, 1);
  model_v_proposed.resize(nparam, 1);

  LeastSquaresIterator::copy(reference, model_0, model_mask);

  int iterations = 0;

  do {
      
    //
    // First copy parameters to model vector
    //
    LeastSquaresIterator::copy(model, model_v, model_mask);

    int m = iterations % 2;
    bool valid = false;
    
    do {
      if (!step[m]->ComputeStepJoint(epsilon[m],
				     Cd_love,
				     Cd_rayleigh,
				     Cm,
				     residuals_love,
				     residuals_rayleigh,
				     G_love,
				     G_rayleigh,
				     dLdp_love,
				     model_mask,
				     model_v,
				     model_0,
				     model_v_proposed)) {
	fprintf(stderr, "error: failed to compute step\n");
	return false;
      }

      valid = LeastSquaresIterator::validate(model_v_proposed, model_mask,
					     PRIOR_MIN,
					     PRIOR_MAX);
      if (!valid) {
	if (epsilon[m] < EPSILON_MIN) {
	  printf("%4d: Prior violation exit\n", iterations);
	  break;
	}
	
	epsilon[m] *= 0.5;
      }
    } while (!valid);

    if (valid) {
      
      LeastSquaresIterator::copy(model_v_proposed, model);
      
      //
      // Recompute Likelihood
      //
      last_like = like;

      if (skip <= 1) {
	like_love = likelihood_love(data_love,
				    model,
				    reference,
				    damping,
				    posterior,
				    mesh,
				    love,
				    dkdp_love,
				    dUdp_love,
				    dLdp_love,
				    G_love,
				    residuals_love,
				    Cd_love,
				    threshold,
				    order,
				    highorder,
				    boundaryorder,
				    scale,
				    frequency_thin);
	
	like_rayleigh = likelihood_rayleigh(data_rayleigh,
					    model,
					    reference,
					    damping,
					    posterior,
					    mesh,
					    rayleigh,
					    dkdp_rayleigh,
					    dUdp_rayleigh,
					    dLdp_rayleigh,
					    G_rayleigh,
					    residuals_rayleigh,
					    Cd_rayleigh,
					    threshold,
					    order,
					    highorder,
					    boundaryorder,
					    scale,
					    frequency_thin);
      } else {
	like_love = likelihood_love_spline(data_love,
					   model,
					   reference,
					   damping,
					   false,
					   mesh,
					   love,
					   dkdp_love,
					   dUdp_love,
					   dLdp_love,
					   G_love,
					   Gk_love,
					   GU_love,
					   residuals_love,
					   Cd_love,
					   threshold,
					   order,
					   highorder,
					   boundaryorder,
					   scale,
					   skip);
	
	like_rayleigh = likelihood_rayleigh_spline(data_rayleigh,
						   model,
						   reference,
						   damping,
						   false,
						   mesh,
						   rayleigh,
						   dkdp_rayleigh,
						   dUdp_rayleigh,
						   dLdp_rayleigh,
						   G_rayleigh,
						   Gk_rayleigh,
						   GU_rayleigh,
						   residuals_rayleigh,
						   Cd_rayleigh,
						   threshold,
						   order,
						   highorder,
						   boundaryorder,
						   scale,
						   skip);
      }	
      
      like = like_love + like_rayleigh;
      
      if (like > last_like) {
	
	if (epsilon[m] < EPSILON_MIN) {
	  printf("%4d: Exiting\n", iterations);
	  break;
	}
	
	//
	// Back track and recompute (a little inefficient here)
	//
	printf("%4d: Backtracking %16.9e %16.9e\n", iterations, like, last_like);
	epsilon[m] *= 0.5;

	//
	// Restore previous model
	//
	LeastSquaresIterator::copy(model_v, model);

	//
	// Restore previous gradients etc.
	//
	residuals_love = old_residuals_love;
	residuals_rayleigh = old_residuals_rayleigh;
	G_love = old_G_love;
	G_rayleigh = old_G_rayleigh;
	dLdp_love = old_dLdp_love;
	
	like = last_like;

      } else {
	
	//
	// Accepted new model, store current gradients etc
	//
	old_residuals_love = residuals_love;
	old_residuals_rayleigh = residuals_rayleigh;
	old_G_love = G_love;
	old_G_rayleigh = G_rayleigh;
	
	for (size_t i = 0; i < nparam; i ++) {
	  dLdp_love(i, 0) += dLdp_rayleigh(i, 0);
	}
	
	old_dLdp_love = dLdp_love;

	printf("%4d: %16.9e %16.9e\n", iterations, like, epsilon[m]);
	
	iterations ++;
	
      }
    } else {
      break;
    }
      
  } while (iterations < maxiterations);

  return true;
}

