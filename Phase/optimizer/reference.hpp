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
#ifndef reference_hpp
#define reference_hpp

#include "common.hpp"

class ReferenceModel
{
public:
  
  ReferenceModel()
  {
  }

  bool load_model(const char *filename)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      return false;
    }

    char cell[1024];
    char hs[1024];

    if (fscanf(fp, "%s\n%s\n", cell, hs) != 2) {
      fprintf(stderr, "error: failed to read cell and halfspace types\n");
      return false;
    }
    if (strcmp(cell, "AnisotropicRhoVsXiVpVs") != 0) {
      fprintf(stderr, "error: cell type invalid\n");
      return false;
    }
    if (strcmp(hs, "AnisotropicRhoVsXiVpVsHalfspace") != 0) {
      fprintf(stderr, "error: halfspace type invalid\n");
      return false;
    }
    
    int maxorder;
    if (fscanf(fp, "%d", &maxorder) != 1) {
      fprintf(stderr, "error: failed to read max order\n");
      return false;
    }

    if (!model.read(fp)) {
      fprintf(stderr, "error: failed to parse model\n");
      return false;
    }

    fclose(fp);
    
    
    reference = model;

    return true;
  }
  
  bool load(const char *filename)
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
