#include <exception>

#include <stdio.h>
#include <getopt.h>

#include "spec1d/isotropicvs.hpp"
#include "spec1d/isotropicvshalfspace.hpp"

#include "spec1d/empiricalmodel.hpp"
#include "spec1d/model.hpp"
#include "spec1d/modelloader.hpp"
#include "spec1d/lovematrices.hpp"
#include "spec1d/lobattoprojection.hpp"

constexpr int MAXORDER = 20;
constexpr int BOUNDARYORDER = MAXORDER;

typedef BrocherEmpiricalModel<double> empiricalmodel_t;
typedef IsotropicVs<double, empiricalmodel_t> node_t;
typedef IsotropicVsHalfspace<double, empiricalmodel_t, MAXORDER> boundary_t;
typedef Cell<double, node_t, MAXORDER> cell_t;

typedef Model<double, node_t, boundary_t, MAXORDER> model_t;
typedef Mesh<double, MAXORDER> mesh_t;
typedef LoveMatrices<double, MAXORDER> lovesolver_t;

static char short_options[] = "i:o:f:F:S:H:p:b:t:P:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},

  {"frequency-min", required_argument, 0, 'f'},
  {"frequency-max", required_argument, 0, 'F'},
  {"frequency-samples", required_argument, 0, 'S'},
  {"frequency-hertz", required_argument, 0, 'H'},
  
  {"scale", required_argument, 0, 's'},
  {"order", required_argument, 0, 'p'},
  {"boundaryorder", required_argument, 0, 'b'},
  {"threshold", required_argument, 0, 't'},
  {"high-order", required_argument, 0, 'P'},

  {"help", no_argument, 0, 'h'},
  
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  char *input_file;
  char *output_file;

  double fmin;
  double fmax;
  double hertz;
  int samples;
  
  double threshold;
  int order;
  int highorder;
  int boundaryorder;
  double scale;

  

  //
  // Defaults
  //
  input_file = nullptr;
  output_file = nullptr;

  fmin = 0.0;
  fmax = 0.5;

  threshold = 0.0;
  order = 5;
  highorder = 5;
  boundaryorder = 5;
  scale = 1.0e-4;

  fmin = 1.0/40.0;
  fmax = 1.0/2.0;

  samples = 8192;
  hertz = 2.0;

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
      input_file = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'f':
      fmin = atof(optarg);
      if (fmin <= 0.0) {
        fprintf(stderr, "error: fmin must be positive\n");
        return -1;
      }
      break;

    case 'F':
      fmax = atof(optarg);
      if (fmax <= 0.0) {
        fprintf(stderr, "error: fmax must be positive\n");
        return -1;
      }
      break;

    case 'S':
      samples = atoi(optarg);
      if (samples <= 0 || samples % 2 != 0) {
	fprintf(stderr, "error: samples must be greater than zero and even\n");
	return -1;
      }
      break;

    case 'H':
      hertz = atof(optarg);
      if (hertz <= 0.0) {
	fprintf(stderr, "error: hertz must be greater than 0\n");
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

    default:
      fprintf(stderr, "unknown option %c\n", c);
    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (input_file == nullptr) {
    fprintf(stderr, "error: missing input file paramter\n");
    return -1;
  }

  if (samples % 2 != 0) {
    fprintf(stderr, "error: uneven no. samples not handled\n");
    return -1;
  }

  model_t model;
  std::vector<double> model_errors;
  mesh_t mesh;
  lovesolver_t love;
  
  //
  // Load the model
  //
  FILE *fp = fopen(input_file, "r");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to open %s for reading\n", input_file);
    return -1;
  }

  int nlayers;
  if (fscanf(fp, "%d\n", &nlayers) != 1) {
    fprintf(stderr, "error: failed to read no. layers\n");
    return -1;
  }

  for (int i = 0; i < nlayers; i ++) {

    double thicknesskm;
    int order;
    if (fscanf(fp, "%lf %d ", &thicknesskm, &order) != 2) {
      fprintf(stderr, "error: failed to read line\n");
      return -1;
    }

    if (i < (nlayers - 1)) {
      //
      // Add normal layer to model
      //
      if (thicknesskm <= 0.0) {
	fprintf(stderr, "error: invalid layer thickness %f layer %d/%d\n", thicknesskm, i, nlayers);
	return -1;
      }
      
      model.cells.push_back(cell_t());

      model.cells[i].thickness = thicknesskm * 1.0e3;
      model.cells[i].order[0] = order;

      double vs, vsstd;
      for (int j = 0; j <= order; j ++) {
	if (fscanf(fp, "%lf", &vs) != 1) {
	  fprintf(stderr, "error: failed to parse vs\n");
	  return -1;
	}
	
	model.cells[i].nodes[j] = node_t(vs * 1.0e3);
	model_errors.push_back(vsstd * 1.0e3);
      }
      if (fscanf(fp, "%lf", &vsstd) != 1) {
	fprintf(stderr, "error: failed to parse vsstd\n");
	return -1;
      }

    } else {
      //
      // Add halfspace layer to model
      //
      if (thicknesskm != 0.0 || order != 0) {
	fprintf(stderr, "error: halfspace layer thickness not zero: %10.6f\n", thicknesskm);
	return -1;
      }

      double vs, vsstd;
      if (fscanf(fp, "%lf %lf", &vs, &vsstd) != 2) {
	fprintf(stderr, "error: failed to parse vs, vsstd\n");
	return -1;
      }

      model.boundary = boundary_t(vs * 1.0e3);
      model_errors.push_back(vsstd * 1.0e3);
      
    }

  }

  fclose(fp);

  //
  // Compute phase and group velocities and Jacobians and linearized estimates of
  // phase/group velocity uncertainties in reverse
  //
  int fsamples = samples/2 + 1;

  double *f = new double[fsamples];
  double *phase = new double[fsamples];
  double *phasestd = new double[fsamples];
  double *group = new double[fsamples];
  double *groupstd = new double[fsamples];

  
  model.project_gradient(mesh, order);

  love.recompute(mesh, boundaryorder, scale);

  int nparameters = -1;
  
  Spec1DMatrix<double> dkdp;
  Spec1DMatrix<double> dUdp;
  Spec1DMatrix<double> dgxdv;
  Spec1DMatrix<double> dgzdv;
  Spec1DMatrix<double> dGvdp;

  for (int i = (fsamples - 1); i >= 0; i --) {

    printf("%d\n", i);
    
    f[i] = (double)i * hertz/(double)samples;
    phase[i] = 0.0;
    phasestd[i] = 0.0;
    group[i] = 0.0;
    groupstd[i] = 0.0;
    
    if (f[i] >= fmin && f[i] <= fmax) {


      double omega = f[i] * 2.0 * M_PI;

      double normA, normB, normC;

      double k = love.solve_fundamental_gradient_sep(mesh,
						     boundaryorder,
						     omega,
						     dkdp,
						     dUdp,
						     normA,
						     normB,
						     normC);

      if (k == 0.0) {
	fprintf(stderr, "error: failed to compute wave number\n");
	return -1;
      }

      if (nparameters < 0) {
	nparameters = dkdp.rows();
	if (nparameters != (int)model_errors.size()) {
	  fprintf(stderr, "error: unexpected no. parameters: %d != %d\n",
		  nparameters, (int)model_errors.size());
	  return -1;
	}
      }
      
      //
      // Compute phase velocity and linearized estimates of errors
      //
      phase[i] = omega/fabs(k);
      for (int j = 0; j < nparameters; j ++) {
	double dcdp = -dkdp(j, 0) * omega/(k*k);
	
	phasestd[i] += (model_errors[j]*model_errors[j]) * (dcdp*dcdp);
      }
      phasestd[i] = sqrt(phasestd[i]);

      //
      // Compute group velocity and linerized estimates of errors
      //
      group[i] = (normB*k)/(omega*normA);
      if (k < 0.0) {
	group[i] = -group[i];
      }
      for (int j = 0; j < nparameters; j ++) {
	groupstd[i] += (model_errors[j]*model_errors[j]) * (dUdp(j, 0)*dUdp(j, 0));
      }
      groupstd[i] = sqrt(groupstd[i]);

      //
      // Update Laguerre scale
      //
      double vs2 = mesh.boundary.L/mesh.boundary.rho;
      
      double disc1 = k*k - omega*omega/vs2;
      
      if (isnormal(disc1) && disc1 > 0.0) {
	
	double newscale1 = sqrt(disc1);
	
	love.recompute(mesh, boundaryorder, newscale1);
      }	
      
    }

  }


  //
  // Save output file
  //

  fp = fopen(output_file, "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create %s\n", output_file);
    return -1;
  }
  
  fprintf(fp, "%d\n", fsamples);

  for (int i = 0; i < fsamples; i ++) {
    fprintf(fp, "%15.9f %16.9e %16.9e %16.9e %16.9e\n",
	    f[i], phase[i], phasestd[i], group[i], groupstd[i]);
  }

  fclose(fp);

  return 0;
}

static void usage(const char *pname)
{
}
