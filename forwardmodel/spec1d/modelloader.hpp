#pragma once
#ifndef modelloader_hpp
#define modelloader_hpp

#include <stdio.h>
#include <string.h>

#include "model.hpp"
#include "logging.hpp"

#include "empiricalmodel.hpp"

#include "isotropicvs.hpp"
#include "isotropicrhovs.hpp"
#include "isotropicrhovpvs.hpp"
#include "anisotropicrhovsxivpvs.hpp"

#include "fixedboundary.hpp"
#include "isotropicvshalfspace.hpp"
#include "isotropicrhovshalfspace.hpp"
#include "isotropicrhovpvshalfspace.hpp"
#include "anisotropicrhovsxivpvshalfspace.hpp"

#include "regression.hpp"

template
<
  typename real,
  size_t maxorder
>
ModelInterface<real, maxorder>* ModelLoader(const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("Failed to open file");
    return nullptr;
  }

  int order;
  if (fscanf(fp, "%d\n", &order) != 1) {
    ERROR("Failed to read header");
    return nullptr;
  }

  if (order > (int)maxorder) {
    ERROR("Order out of range");
    return nullptr;
  }
  
  char parameterset[1024];
  char boundary[1024];

  if (fscanf(fp, "%s\n%s\n", parameterset, boundary) != 2) {
    ERROR("Failed to read types");
    return nullptr;
  }

  ModelInterface<real, maxorder> *p = nullptr;
  
  if (strcmp(parameterset, "IsotropicRhoVpVs") == 0) {
    if (strcmp(boundary, "FixedBoundary") == 0) {

      p = new Model<real, IsotropicRhoVpVs<real>, FixedBoundary<real, maxorder>, maxorder>();

    } else if (strcmp(boundary, "IsotropicRhoVpVsHalfspace") == 0) {

      p = new Model<real, IsotropicRhoVpVs<real>, IsotropicRhoVpVsHalfspace<real, maxorder>, maxorder>();

    } else {
      ERROR("Unknown/incompatible boundary: %s", boundary);
      return nullptr;
    }
    
  } else if (strcmp(parameterset, "IsotropicVs") == 0) {

    if (strcmp(boundary, "FixedBoundary") == 0) {

      p = new Model<real,
		    IsotropicVs<real, BrocherEmpiricalModel<real>>,
		    FixedBoundary<real, maxorder>,
		    maxorder>();

    } else if (strcmp(boundary, "IsotropicVsHalfspace") == 0) {

      p = new Model<real,
		    IsotropicVs<real, BrocherEmpiricalModel<real>>,
		    IsotropicVsHalfspace<real, BrocherEmpiricalModel<real>, maxorder>,
		    maxorder>();

    } else {
      ERROR("Unknown/incompatible boundary: %s", boundary);
      return nullptr;
    }

  } else if (strcmp(parameterset, "AnisotropicRhoVsXiVpVs") == 0) {

    if (strcmp(boundary, "FixedBoundary") == 0) {

      p = new Model<real,
		    AnisotropicRhoVsXiVpVs<real>,
		    FixedBoundary<real, maxorder>,
		    maxorder>();
      
    } else if (strcmp(boundary, "AnisotropicRhoVsXiVpVsHalfspace") == 0) {

      p = new Model<real,
		    AnisotropicRhoVsXiVpVs<real>,
		    AnisotropicRhoVsXiVpVsHalfspace<real, maxorder>,
		    maxorder>();

    } else {
      ERROR("Unknown/incompatible boundary: %s", boundary);
      return nullptr;
    }
    
  } else if (strcmp(parameterset, "Regression") == 0) {

    if (strcmp(boundary, "FixedBoundary") == 0) {

      p = new Model<real, Regression<real>, FixedBoundary<real, maxorder>, maxorder>();

    } else {
      ERROR("Unknown/incompatible boundary: %s", boundary);
      return nullptr;
    }
      
  } else {

    ERROR("Unknown parameter set: %s", parameterset);
    return nullptr;
  }

  if (p == nullptr) {
    ERROR("Null model");
    return nullptr;
  }

  if (!p->read(fp)) {
    ERROR("Failed to read model");
    return nullptr;
  }

  fclose(fp);
  return p;
}

#endif // modelloaded_hpp
  
