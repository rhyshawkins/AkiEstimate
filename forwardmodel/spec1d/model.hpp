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
#ifndef model_hpp
#define model_hpp

#include <vector>

#include "modelinterface.hpp"
#include "cell.hpp"

#include "mesh.hpp"

#include "encodedecode.hpp"

template
<
  typename real,
  typename parameterset,
  typename boundarycondition,
  size_t maxorder
>
class Model : public ModelInterface<real, maxorder> {
public:

  Model()
  {
  }

  virtual ~Model()
  {
  }
  
  bool read(FILE *fp)
  {
    int fcells;
    if (fscanf(fp, "%d\n", &fcells) != 1) {
      ERROR("Failed to read no. cells");
      return false;
    }

    if (fcells < 0) {
      ERROR("No. cells out of range");
      return false;
    }

    cells.clear();
    cells.resize(fcells);
    
    for (int i = 0; i < fcells; i ++) {
      if (!cells[i].read(fp)) {
	ERROR("Failed to read cell");
	return false;
      }
    }

    if (!boundary.read(fp)) {
      ERROR("Failed to read boundary");
      return false;
    }

    return true;
  }

  virtual void project(Mesh<real, maxorder> &mesh, size_t order) const
  {
    mesh.cells.clear();

    real offset = 0.0;
    size_t i = 0;
    
    for (auto &c : cells) {
      c.project(i, offset, mesh, order);
      offset += c.thickness;
      i ++;
    }
      
    boundary.project(i, offset, mesh);
  }

  virtual void project_gradient(Mesh<real, maxorder> &mesh, size_t order) const
  {
    mesh.cells.clear();
    mesh.cell_parameter_offsets.clear();
    mesh.cell_parameter_offsets.push_back(0);

    mesh.cell_thickness_reference.clear();
    
    real offset = 0.0;

    size_t i = 0;
    for (auto &c : cells) {
      c.project_gradient(i, offset, mesh, order);
      offset += c.thickness;
      i ++;
    }
      
    boundary.project_gradient(i, offset, mesh);
  }

  virtual void project_threshold(Mesh<real, maxorder> &mesh,
				 double threshold,
				 size_t high_order,
				 double high_refinement,
				 size_t low_order)
  {
    mesh.cells.clear();

    real offset = 0.0;

    size_t i = 0;
    for (auto &c : cells) {

      if (threshold <= 0.0 || offset > threshold) {
	c.project(i, offset, mesh, low_order);
      } else {
	c.project_threshold(i, offset, mesh, threshold, high_refinement, high_order, low_order);
      }
      
      offset += c.thickness;
      i ++;
    }
      
    boundary.project(i, offset, mesh);
  }
				 
  virtual void project_threshold_gradient(Mesh<real, maxorder> &mesh,
					  double threshold,
					  size_t high_order,
					  double high_refinement,
					  size_t low_order)
  {
    mesh.cells.clear();
    mesh.cell_parameter_offsets.clear();
    mesh.cell_parameter_offsets.push_back(0);

    real offset = 0.0;

    size_t i = 0;
    for (auto &c : cells) {

      if (threshold <= 0.0 || offset > threshold) {
	c.project_gradient(i, offset, mesh, low_order);
      } else {
	c.project_threshold_gradient(i, offset, mesh, threshold, high_refinement, high_order, low_order);
      }
      
      offset += c.thickness;
      i ++;
    }
      
    boundary.project_gradient(i, offset, mesh);
  }

  virtual real halfspace_depth() const
  {
    real depth = 0.0;

    for (auto &c: cells) {
      depth += c.thickness;
    }

    return depth;
  }
  
  virtual real compute_node_depth(const Mesh<real, maxorder> &mesh, int cell, int order, int node) const
  {
    real offset = 0.0;

    if (cell < 0 || cell >= (int)cells.size()) {
      FATAL("Invalid cell: %d (%d)", cell, (int)cells.size());
    }

    if (order < 0 || order > (int)maxorder) {
      FATAL("Invalid order %d (%d)", order, (int)maxorder);
    }

    if (node < 0 || node > order) {
      FATAL("Invalid node");
    }
    
    for (int i = 0; i < cell; i ++) {
      offset += cells[i].thickness;
    }

    real xi = mesh.quadrature[order]->nodes[node];

    return offset + cells[cell].thickness * (xi + 1.0)/2.0;
  }

  virtual bool interpolate(const Mesh<real, maxorder> &mesh, real z, real *parameters) const
  {
    real offset = 0.0;
      
    for (auto &c : cells) {

      if (offset + c.thickness >= z) {

	real xi = 2.0 * (z - offset)/c.thickness - 1.0;
	
	for (int j = 0; j < (int)parameterset::NPARAMETERS; j ++) {

	  parameters[j] = 0.0;

	  for (int i = 0; i <= (int)c.order[j]; i ++) {

	    parameters[j] += mesh.quadrature[c.order[j]]->cardinal[i].value(xi) * c.nodes[i][j];

	  }

	}

	return true;
      }

      offset += c.thickness;
    }

    for (int j = 0; j < (int)parameterset::NPARAMETERS; j ++) {
      parameters[j] = boundary.parameters[j];
    }
    return true;
  }

  virtual size_t parameter_count() const
  {
    size_t count = boundary.nparameters();
    
    for (auto &c : cells) {
      for (int j = 0; j < (int)parameterset::NPARAMETERS; j ++) {

	count += (c.order[j] + 1);
	
      }
    }

    return count;
  }
  
  virtual bool interpolate_gradient(const Mesh<real, maxorder> &mesh,
				    real z,
				    real *parameters,
				    Spec1DMatrix<real> &dvdp) const
  {
    real offset = 0.0;
    size_t poffset = 0;
      
    for (auto &c : cells) {

      if (offset + c.thickness >= z) {

	real xi = 2.0 * (z - offset)/c.thickness - 1.0;
	
	for (int j = 0; j < (int)parameterset::NPARAMETERS; j ++) {

	  parameters[j] = 0.0;

	  for (int i = 0; i <= (int)c.order[j]; i ++) {

	    real w = mesh.quadrature[c.order[j]]->cardinal[i].value(xi);
	    
	    parameters[j] += w * c.nodes[i][j];

	    dvdp(poffset + i, 0) = w;

	  }

	  poffset += (c.order[j] + 1);

	}

	return true;

      } else {
	for (int j = 0; j < (int)parameterset::NPARAMETERS; j ++) {
	  poffset += (c.order[j] + 1);
	}
      }

      offset += c.thickness;
    }

    for (int j = 0; j < (int)parameterset::NPARAMETERS; j ++) {
      parameters[j] = boundary.parameters[j];
    }
    
    return true;
  }

  virtual void project_with_refinement(Mesh<real, maxorder> &mesh,
				       double min_basement,
				       double max_cell_thickness,
				       size_t order) const
  {
    mesh.cells.clear();

    real offset = 0.0;
    real top_offset = 0.0;
    size_t i = 0;
    for (auto &c : cells) {
      top_offset = offset;
      c.project_with_refinement(i, offset, mesh, max_cell_thickness, order);
      offset += c.thickness;
      i ++;
    }

    MeshParameter<real> lower_boundary = cells[cells.size() - 1].lower_boundary(top_offset);
    
    boundary.project(i,
		     offset,
		     min_basement,
		     max_cell_thickness,
		     lower_boundary,
		     mesh);
  }

  void print(FILE *fp)
  {
    for (auto &c: cells) {
      c.print(fp);
    }

    boundary.print(fp);
  }

  bool save(const char *filename) const
  {
    FILE *fp = fopen(filename, "w");

    if (fp == NULL) {
      ERROR("Failed to create file");
      return false;
    }

    fprintf(fp, "%s\n", parameterset::NAME());
    fprintf(fp, "%s\n", boundarycondition::NAME());
    
    fprintf(fp, "%d %d\n", (int)maxorder, (int)cells.size());

    for (auto &c: cells) {
      if (!c.save(fp)) {
	ERROR("Failed to save cell\n");
	return false;
      }
    }

    if (!boundary.save(fp)) {
      ERROR("Failed to save boundary\n");
      return false;
    }

    fclose(fp);
    return true;
  }

  int encode_size()
  {
    int size = sizeof(int);

    for (auto &c : cells) {
      size += c.encode_size();
    }

    size += boundary.encode_size();
    return size;
  }
  
  int encode(char *buffer, int &offset, int buffer_size)
  {
    int s = 0;
    int e;
    int c = cells.size();

    e = ::encode<int>(c, buffer, offset, buffer_size);
    if (e < 0) {
      ERROR("Failed to encode no. cells");
      return -1;
    }
    s += e;

    for (auto &c : cells) {
      e = c.encode(buffer, offset, buffer_size);
      if (e < 0) {
	ERROR("Failed to encode cell %d", (int)cells.size());
	return -1;
      }
      s += e;
    }

    e = boundary.encode(buffer, offset, buffer_size);
    if (e < 0) {
      ERROR("Failed to encode boundary");
      return -1;
    }
    s += e;

    return s;
  }

  int decode(char *buffer, int &offset, int buffer_size)
  {
    int c;

    if (::decode<int>(c, buffer, offset, buffer_size) < 0) {
      ERROR("Failed to decode number of cells");
      return -1;
    }

    if (c <= 0) {
      ERROR("Invalid no. of cells: %d", c);
      return -1;
    }
    
    cells.resize(c);

    for (auto &c : cells) {
      
      if (c.decode(buffer, offset, buffer_size) < 0) {
	ERROR("Failed to decode cell");
	return -1;
      }
    }

    return boundary.decode(buffer, offset, buffer_size);
  }

  
  bool load(const char *filename)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      ERROR("Failed to open file");
      return false;
    }

    int fmaxorder, fcells;
    if (fscanf(fp, "%d %d\n", &fmaxorder, &fcells) != 2) {
      ERROR("Failed to read header");
      return false;
    }

    if (fmaxorder < 0 || fmaxorder > (int)maxorder) {
      ERROR("Max order out of range");
      return false;
    }

    if (fcells < 0) {
      ERROR("No. cells out of range");
      return false;
    }

    cells.clear();
    cells.resize(fcells);
    
    for (int i = 0; i < fcells; i ++) {
      if (!cells[i].load(fp)) {
	ERROR("Failed to read cell");
	return false;
      }
    }

    if (!boundary.load(fp)) {
      ERROR("Failed to read boundary");
      return false;
    }

    fclose(fp);
    return true;
  }

  void clone_with_split(const Model &model, size_t cell, double split)
  {
    cells.clear();
    
    for (size_t i = 0; i < model.cells.size(); i ++) {

      if (i == cell) {

	cells.push_back(Cell<real,parameterset,maxorder>());
	cells.push_back(Cell<real,parameterset,maxorder>());

	cells[i].thickness = (1.0 - split) * model.cells[i].thickness;
	cells[i + 1].thickness = split * model.cells[i].thickness;

      } else {
	cells.push_back(model.cells[i]);
      }
    }

    boundary = model.boundary;
  }

  void clone_with_merge(const Model &model, size_t cell, double &split)
  {
    cells.clear();

    for (size_t i = 0; i < model.cells.size(); i ++) {

      if (i == cell) {

	if (i == (model.cells.size() - 1)) {
	  FATAL("i out of range");
	}
	
	real thicka = model.cells[i].thickness;
	real thickb = model.cells[i + 1].thickness;

	cells.push_back(Cell<real,parameterset,maxorder>());

	cells[i].thickness = thicka + thickb;
	split = thickb/thicka + thickb;

      } else if (i == (cell + 1)) {
	// Do nothing, this cell is removed
      } else {
	cells.push_back(model.cells[i]);
      }
    }

    boundary = model.boundary;
  }

  real totalthickness() const
  {
    real t = 0;
    for (auto &c: cells) {
      t += c.thickness;
    }
    return t;
  }
  
  std::vector<Cell<real, parameterset, maxorder>> cells;
  boundarycondition boundary;
};

#endif // model_hpp
