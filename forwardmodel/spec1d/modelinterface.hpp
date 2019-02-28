#pragma once
#ifndef modelinterface_hpp
#define modelinterface_hpp

#include "mesh.hpp"

template
<
  typename real,
  size_t maxorder
>
class ModelInterface {
public:

  ModelInterface()
  {
  }

  virtual ~ModelInterface()
  {
  }
  
  virtual bool read(FILE *fp) = 0;
  
  virtual void project(Mesh<real, maxorder> &mesh, size_t order) const = 0;
  virtual void project_with_refinement(Mesh<real, maxorder> &mesh,
				       double min_basement,
				       double max_cell_thickness,
				       size_t order) const = 0;
  virtual void project_threshold(Mesh<real, maxorder> &mesh,
                                 double threshold,
                                 size_t high_order,
                                 double high_refinement,
                                 size_t low_order) = 0;

  virtual bool interpolate(const Mesh<real, maxorder> &mesh, real z, real *parameters) const = 0;

  virtual size_t parameter_count() const = 0;
  
  virtual bool interpolate_gradient(const Mesh<real, maxorder> &mesh,
				    real z,
				    real *parameters,
				    Spec1DMatrix<real> &dvdp) const = 0;
  
};

#endif // modelinterface_hpp
