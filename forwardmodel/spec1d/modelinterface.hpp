//
//    Spec1D : A spectral element code for surface wave dispersion of Love
//    and Rayleigh waves. See
//
//      R Hawkins, "A spectral element method for surface wave dispersion and adjoints",
//      Geophysical Journal International, 2018, 215:1, 267 - 302
//      https://doi.org/10.1093/gji/ggy277
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
