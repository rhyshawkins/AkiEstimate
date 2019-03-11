#pragma once
#ifndef reference_hpp
#define reference_hpp

#include "common.hpp"

class ReferenceModel
{
public:
  
  ReferenceModel()
  {
  }

  bool load(const char *filename, bool promote, size_t promote_order)
  {

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to open %s for reading\n", filename);
      return -1;
    }
    
    int nlayers;
    if (fscanf(fp, "%d\n", &nlayers) != 1) {
      fprintf(stderr, "error: failed to read no. layers\n");
      return false;
    }
    
    for (int i = 0; i < nlayers; i ++) {
      
      double thicknesskm;
      int corder;
      double vs;
      double vsstd;
      if (fscanf(fp, "%lf %d", &thicknesskm, &corder) != 2) {
	fprintf(stderr, "error: failed to read line\n");
	return false;
      }

      if (i < (nlayers - 1)) {
	  
	//
	// Add normal layer to model
	//
	if (thicknesskm <= 0.0) {
	  fprintf(stderr, "error: invalid layer thickness %f layer %d/%d\n", thicknesskm, i, nlayers);
	  return false;
	}
	
	if (corder == 0 && promote) {

	  if (fscanf(fp, "%lf", &vs) != 1) {
	    fprintf(stderr, "error: failed to read vs\n");
	    return -1;
	  }
	  
	  double rho;
	  double vp;
	  vs *= 1.0e3;
	  BrocherEmpiricalModel<double>::compute(vs, rho, vp);
	  double vpvs = vp/vs;

	  model.cells.push_back(cell_t());
	  model.cells[i].thickness = thicknesskm * 1.0e3;
	  model.cells[i].order[0] = promote_order;
	  model.cells[i].order[1] = promote_order;
	  model.cells[i].order[2] = promote_order;
	  model.cells[i].order[3] = promote_order;

	  for (size_t j = 0; j <= promote_order; j ++) {
	    model.cells[i].nodes[j] = cell_parameter_t(rho, vs, 1.0, vpvs);
	  }

	  if (fscanf(fp, "%lf", &vsstd) != 1) {
	    fprintf(stderr, "error: failed to read vs error\n");
	    return -1;
	  }
	  
	} else {
	  
	  model.cells.push_back(cell_t());
	  model.cells[i].thickness = thicknesskm * 1.0e3;
	  model.cells[i].order[0] = corder;
	  model.cells[i].order[1] = corder;
	  model.cells[i].order[2] = corder;
	  model.cells[i].order[3] = corder;
	  
	  for (int j = 0; j <= corder; j ++) {
	    
	    if (fscanf(fp, "%lf", &vs) != 1) {
	      fprintf(stderr, "error: failed to read vs\n");
	      return -1;
	    }
	    
	    double rho;
	    double vp;
	    vs *= 1.0e3;
	    BrocherEmpiricalModel<double>::compute(vs, rho, vp);
	    double vpvs = vp/vs;
	    
	    model.cells[i].nodes[j] = cell_parameter_t(rho, vs, 1.0, vpvs);
	  }
	  
	  if (fscanf(fp, "%lf", &vsstd) != 1) {
	    fprintf(stderr, "error: failed to read vs error\n");
	    return -1;
	  }
	}
	
      } else {
	if (fscanf(fp, "%lf %lf", &vs, &vsstd) != 2) {
	  fprintf(stderr, "error: failed to read vs\n");
	  return -1;
	}
	
	//
	// Add halfspace layer to model
	//
	double rho;
	double vp;
	vs *= 1.0e3;
	BrocherEmpiricalModel<double>::compute(vs, rho, vp);
	double vpvs = vp/vs;
	
	if (thicknesskm != 0.0) {
	  fprintf(stderr, "error: halfspace layer thickness not zero: %10.6f\n", thicknesskm);
	  return -1;
	}
	
	model.boundary = boundary_t(rho, vs, 1.0, vpvs);
	
      }
    }
      
    fclose(fp);

    reference = model;
    return true;
    
  }

  model_t model;
  model_t reference;
};

#endif // reference_hpp
