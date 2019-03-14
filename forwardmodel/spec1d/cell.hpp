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
#ifndef cell_hpp
#define cell_hpp

#include <array>

#include "mesh.hpp"
#include "lobattoprojection.hpp"

#include "encodedecode.hpp"

template
<
  typename real,
  typename parameterset,
  size_t maxorder
>
class Cell {
public:

  Cell()
  {
    for (size_t i = 0; i < parameterset::NPARAMETERS; i ++) {
      order[i] = 0;
    }
  }

  ~Cell()
  {
  }

  size_t nparameters() const
  {
    return parameterset::NPARAMETERS;
  }

  void node_project(Mesh<real, maxorder> &mesh, size_t mesh_order) const
  {
    //
    // Project current order(s) for each parameter to specified order (nodes -> projected nodes)
    //
    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      for (size_t j = 0; j <= mesh_order; j ++) {
	projected_nodes[j][k] = 0.0;

	for (size_t i = 0; i <= order[k]; i ++) {

	  projected_nodes[j][k] += nodes[i][k] * mesh.projection[mesh_order]->weight[order[k]][i][j];

	}
      }
    }
  }
  
  void node_project_gradient(real offset,
			     real thickness,
			     Mesh<real, maxorder> &mesh,
			     Spec1DMatrix<real> &jacobian,
			     size_t mesh_order) const
  {
    //
    // Project current order(s) for each parameter to specified order (nodes -> projected nodes)
    //
    size_t poffset = 0;
    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      for (size_t j = 0; j <= mesh_order; j ++) {
	projected_nodes[j][k] = 0.0;

	real depth = offset + ((mesh.quadrature[maxorder]->nodes[k] + 1.0)/2.0) * thickness;
	
	for (size_t i = 0; i <= order[k]; i ++) {

	  double w = mesh.projection[mesh_order]->weight[order[k]][i][j];
	  
	  projected_nodes[j][k] += nodes[i][k] * w;

	  //
	  // Col offset is 6 * j
	  // Row offset is dynamic = \sum_k (order_i + 1) == poffset
	  //
	  jacobian(poffset + i, 6*j + 0) += w * nodes[i].drho(k, depth);
	  jacobian(poffset + i, 6*j + 1) += w * nodes[i].dA(k, depth);
	  jacobian(poffset + i, 6*j + 2) += w * nodes[i].dC(k, depth);
	  jacobian(poffset + i, 6*j + 3) += w * nodes[i].dF(k, depth);
	  jacobian(poffset + i, 6*j + 4) += w * nodes[i].dL(k, depth);
	  jacobian(poffset + i, 6*j + 5) += w * nodes[i].dN(k, depth);
	}
      }

      poffset += (order[k] + 1);
    }
  }

  void project(size_t index, real offset, Mesh<real, maxorder> &mesh, size_t mesh_order) const
  {
    MeshCell<real, maxorder> mcell(index, thickness);

    mcell.order = mesh_order;
    
    //
    // Project cell from current order to mesh order: nodes -> projected_nodes
    //
    node_project(mesh, mesh_order);

    //
    // Project cell parameters to 6 TI parameters
    //
    for (size_t k = 0; k <= mesh_order; k ++) {

      real depth = offset + ((mesh.quadrature[maxorder]->nodes[k] + 1.0)/2.0) * thickness;
      
      mcell.nodes[k].rho += projected_nodes[k].rho(depth);
      mcell.nodes[k].A += projected_nodes[k].A(depth);
      mcell.nodes[k].C += projected_nodes[k].C(depth);
      mcell.nodes[k].F += projected_nodes[k].F(depth);
      mcell.nodes[k].L += projected_nodes[k].L(depth);
      mcell.nodes[k].N += projected_nodes[k].N(depth);
      
    }
    
    mesh.cells.push_back(mcell);
  }

  void project_gradient(size_t index, real offset, Mesh<real, maxorder> &mesh, size_t mesh_order) const
  {
    MeshCell<real, maxorder> mcell(index, thickness);

    size_t nparameters = 0;
    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {
      nparameters += (order[k] + 1);
    }

    
    mesh.cell_parameter_offsets.push_back(mesh.cell_parameter_offsets[mesh.cell_parameter_offsets.size() - 1] +
					  nparameters);

    if (mesh.cell_thickness_reference.size() == 0) {
      mesh.cell_thickness_reference.push_back(0);
    } else { 
      mesh.cell_thickness_reference.push_back(mesh.cell_thickness_reference[mesh.cell_thickness_reference.size() - 1]);
    }

    mcell.thickness_jacobian = 1.0;
    
    mcell.order = mesh_order;
    mcell.jacobian.resize(nparameters, (mesh_order + 1) * 6);
    mcell.jacobian.setZero();

    node_project_gradient(offset, thickness, mesh, mcell.jacobian, mesh_order);

    //
    // Project cell parameters to 6 TI parameters
    //
    for (size_t k = 0; k <= mesh_order; k ++) {

      real depth = offset + ((mesh.quadrature[maxorder]->nodes[k] + 1.0)/2.0) * thickness;
      
      mcell.nodes[k].rho += projected_nodes[k].rho(depth);
      mcell.nodes[k].A += projected_nodes[k].A(depth);
      mcell.nodes[k].C += projected_nodes[k].C(depth);
      mcell.nodes[k].F += projected_nodes[k].F(depth);
      mcell.nodes[k].L += projected_nodes[k].L(depth);
      mcell.nodes[k].N += projected_nodes[k].N(depth);
      
    }

    mesh.cells.push_back(mcell);

  }

  void cell_project(LobattoUpProjector<real, maxorder> &up_projector,
		    LobattoDownProjector<real, maxorder> &down_projector,
		    Cell &cell) const
  {
    for (size_t k = 0; k < cell.nodes[0].parameters.size(); k ++) {

      if (cell.order[k] == order[k]) {

	//
	// Copy
	//
	for (size_t j = 0; j <= cell.order; j ++) {
	  cell.nodes[j].parameters[k] = nodes[j].parameters[k];
	}
	
      } else if (cell.order[k] == (order[k] + 1)) {

	//
	// Up projected
	//
	for (size_t j = 0; j <= cell.order; j ++) {
	  cell.nodes[j].parameters[k] = 0.0;
	}

	for (size_t i = 0; i <= order; i ++) {
	  for (size_t j = 0; j <= cell.order; j ++) {
	    cell.nodes[j].parameters[k] += nodes[i].parameters[k] * up_projector.P[i].weights[i][j];
	  }
	}
	
      } else if ((cell.order[k] + 1) == order[k]) {

	//
	// Down projected
	//
	for (size_t j = 0; j <= cell.order; j ++) {
	  cell.nodes[j].parameters[k] = 0.0;
	}
	for (size_t i = 0; i <= order; i ++) {
	  for (size_t j = 0; j <= cell.order; j ++) {
	    cell.nodes[j].parameters[k] += nodes[i].parameters[k] * down_projector.P[i].weights[i][j];
	  }
	}
      } else {

	FATAL("Only +-1 order transitions are implemented");

      }
    }
	
  }

  void project_threshold(size_t index,
			 real offset,
			 Mesh<real, maxorder> &mesh,
			 double threshold,
			 double high_refinement,
			 size_t high_order,
			 size_t low_order)
  {
    if (offset >= threshold) {
      //
      // Cell completely below threshold
      //
      project(index, offset, mesh, low_order);
      
    } else if ((offset + thickness) <= threshold) {
      //
      // Cell completely above threshold
      //
      project_with_refinement(index,
			      offset,
			      mesh,
			      high_refinement,
			      high_order);

    } else {
      //
      // Cell cut by threshold
      //
      
      if (high_refinement <= 0.0 || high_refinement > thickness) {

	// Split cell at threshold into high and low order
	// node_project(mesh, high_order);
	
	real cell1_thickness = threshold - offset;
	MeshCell<real, maxorder> cell1(index, cell1_thickness);
	cell1.order = high_order;
	
	for (size_t j = 0; j <= high_order; j ++) {
	  //
	  // Determine the xi of the sub cell node in the cell
	  //
	  real zp = (mesh.quadrature[high_order]->nodes[j] + 1.0)/2.0 * cell1_thickness;
	  real xi = 2.0 * zp/thickness - 1.0;

	  parameterset p = interpolate(mesh, xi);

	  real depth = offset + zp;
	  cell1.nodes[j].rho = p.rho(depth);
	  cell1.nodes[j].A = p.A(depth);
	  cell1.nodes[j].C = p.C(depth);
	  cell1.nodes[j].F = p.F(depth);
	  cell1.nodes[j].L = p.L(depth);
	  cell1.nodes[j].N = p.N(depth);
	}
	mesh.cells.push_back(cell1);

	real cell2_thickness = thickness - cell1_thickness;
	offset += cell1_thickness;

	MeshCell<real, maxorder> cell2(index, cell2_thickness);
	cell2.order = low_order;
	
	for (size_t j = 0; j <= low_order; j ++) {
	  //
	  // Determine the xi of the sub cell node in the cell
	  //
	  real zp = cell1_thickness + (mesh.quadrature[low_order]->nodes[j] + 1.0)/2.0 * cell2_thickness;
	  real xi = 2.0 * zp/thickness - 1.0;
	  
	  parameterset p = interpolate(mesh, xi);

	  real depth = offset + zp;
	  cell2.nodes[j].rho = p.rho(depth);
	  cell2.nodes[j].A = p.A(depth);
	  cell2.nodes[j].C = p.C(depth);
	  cell2.nodes[j].F = p.F(depth);
	  cell2.nodes[j].L = p.L(depth);
	  cell2.nodes[j].N = p.N(depth);
	}
	mesh.cells.push_back(cell2);
	
      } else {

	FATAL("Unimplemented");

      }
    }
  }

  void project_threshold_gradient(size_t index,
				  real offset,
				  Mesh<real, maxorder> &mesh,
				  double threshold,
				  double high_refinement,
				  size_t high_order,
				  size_t low_order)
  {
    if (offset > threshold) {
      //
      // Cell completely below threshold
      //
      project_gradient(index, offset, mesh, low_order);
      
    } else if ((offset + thickness) < threshold) {
      //
      // Cell completely above threshold
      //
      project_with_refinement_gradient(index,
				       offset,
				       mesh,
				       high_refinement,
				       high_order);
      
    } else {

      size_t nparameters = 0;
      for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {
	nparameters += (order[k] + 1);
      }
      
      //
      // Cell cut by threshold
      //
      if (high_refinement <= 0.0 || high_refinement > thickness) {

	size_t base_offset = mesh.cell_parameter_offsets[mesh.cell_parameter_offsets.size() - 1];
	size_t parameters_offset = base_offset + nparameters;
	
	mesh.cell_parameter_offsets.push_back(base_offset);
	mesh.cell_parameter_offsets.push_back(parameters_offset);

	real cell1_thickness = threshold - offset;
	MeshCell<real, maxorder> cell1(index, cell1_thickness);
	cell1.order = high_order;
	cell1.jacobian.resize(nparameters, (high_order + 1) * 6);
	cell1.jacobian.setZero();

	// Split cell at threshold into high and low order
	// node_project_gradient(offset, cell1_thickness, mesh, cell1.jacobian, high_order);

	for (size_t j = 0; j <= high_order; j ++) {
	  //
	  // Determine the xi of the sub cell node in the cell
	  //
	  real zp = (mesh.quadrature[high_order]->nodes[j] + 1.0)/2.0 * cell1_thickness;
	  real xi = 2.0 * zp/thickness - 1.0;
	  real depth = offset + zp;
	  
	  parameterset p = interpolate_gradient(mesh, xi, depth, j, cell1.jacobian);
	  cell1.nodes[j].rho = p.rho(depth);
	  cell1.nodes[j].A = p.A(depth);
	  cell1.nodes[j].C = p.C(depth);
	  cell1.nodes[j].F = p.F(depth);
	  cell1.nodes[j].L = p.L(depth);
	  cell1.nodes[j].N = p.N(depth);
	  
	}
	mesh.cells.push_back(cell1);

	real cell2_thickness = thickness - cell1_thickness;
	offset += cell1_thickness;
	MeshCell<real, maxorder> cell2(index, cell2_thickness);
	cell2.order = low_order;
	cell2.jacobian.resize(nparameters, (low_order + 1) * 6);
	cell2.jacobian.setZero();
	
	for (size_t j = 0; j <= low_order; j ++) {
	  //
	  // Determine the xi of the sub cell node in the cell
	  //
	  real zp = cell1_thickness + (mesh.quadrature[low_order]->nodes[j] + 1.0)/2.0 * cell2_thickness;
	  real xi = 2.0 * zp/thickness - 1.0;
	  real depth = offset + zp;
	  
	  parameterset p = interpolate_gradient(mesh, xi, depth, j, cell2.jacobian);

	  cell2.nodes[j].rho = p.rho(depth);
	  cell2.nodes[j].A = p.A(depth);
	  cell2.nodes[j].C = p.C(depth);
	  cell2.nodes[j].F = p.F(depth);
	  cell2.nodes[j].L = p.L(depth);
	  cell2.nodes[j].N = p.N(depth);
	}
	
	mesh.cells.push_back(cell2);
	
      } else {

	FATAL("Unimplemented");

      }
    }
  }

  void project_with_refinement(size_t index,
			       real offset,
			       Mesh<real, maxorder> &mesh,
			       double max_cell_thickness,
			       size_t mesh_order) const
  {
    if (max_cell_thickness <= 0.0 || thickness < max_cell_thickness) {
      project(index, offset, mesh, mesh_order);
    } else {

      node_project(mesh, mesh_order);

      size_t nrefinedcells = ceil(thickness/max_cell_thickness);
      real refinedthickness = thickness/nrefinedcells;

      real refined_offset = 0.0;

      for (size_t i = 0; i < nrefinedcells; i ++) {
      
	MeshCell<real, maxorder> mcell(index, refinedthickness);

	mcell.order = mesh_order;
	
	//
	// Need to interpolate to sub cell.
	//
	for (size_t j = 0; j <= mesh_order; j ++) {

	  //
	  // Determine the xi of the sub cell node in the cell
	  //
	  real zp = refined_offset + (mesh.quadrature[mesh_order]->nodes[j] + 1.0)/2.0 * refinedthickness;
	  real xi = 2.0 * zp/thickness - 1.0;
	  
	  mcell.nodes[j].zero();
	  
	  for (size_t k = 0; k <= mesh_order; k ++) {
	    double l = mesh.quadrature[mesh_order]->cardinal[k].value(xi);

	    real depth = offset + ((xi + 1.0)/2.0 * thickness);
	      
	    mcell.nodes[j].rho += projected_nodes[k].rho(depth) * l;
	    mcell.nodes[j].A += projected_nodes[k].A(depth) * l;
	    mcell.nodes[j].C += projected_nodes[k].C(depth) * l;
	    mcell.nodes[j].F += projected_nodes[k].F(depth) * l;
	    mcell.nodes[j].L += projected_nodes[k].L(depth) * l;
	    mcell.nodes[j].N += projected_nodes[k].N(depth) * l;
	  }
	  
	}
	
	mesh.cells.push_back(mcell);
	refined_offset += refinedthickness;
      }
    }
  }

  void project_with_refinement_gradient(size_t index,
					real offset,
					Mesh<real, maxorder> &mesh,
					double max_cell_thickness,
					size_t mesh_order) const
  {
    if (max_cell_thickness <= 0.0 || thickness < max_cell_thickness) {
      project_gradient(index, offset, mesh, mesh_order);
    } else {

      size_t nrefinedcells = ceil(thickness/max_cell_thickness);
      real refinedthickness = thickness/nrefinedcells;

      real refined_offset = 0.0;

      for (size_t i = 0; i < nrefinedcells; i ++) {
      
	MeshCell<real, maxorder> mcell(index, refinedthickness);

	mcell.order = mesh_order;
	node_project_gradient(refined_offset, refinedthickness, mesh, mcell.jacobian, mesh_order);

	//
	// Need to interpolate to sub cell.
	//
	for (size_t j = 0; j <= mesh_order; j ++) {

	  //
	  // Determine the xi of the sub cell node in the cell
	  //
	  real zp = refined_offset + (mesh.quadrature[mesh_order]->nodes[j] + 1.0)/2.0 * refinedthickness;
	  real xi = 2.0 * zp/thickness - 1.0;
	  
	  mcell.nodes[j].zero();
	  
	  for (size_t k = 0; k <= mesh_order; k ++) {
	    double l = mesh.quadrature[mesh_order]->cardinal[k].value(xi);

	    real depth = offset + ((xi + 1.0)/2.0 * thickness);
	      
	    mcell.nodes[j].rho += projected_nodes[k].rho(depth) * l;
	    mcell.nodes[j].A += projected_nodes[k].A(depth) * l;
	    mcell.nodes[j].C += projected_nodes[k].C(depth) * l;
	    mcell.nodes[j].F += projected_nodes[k].F(depth) * l;
	    mcell.nodes[j].L += projected_nodes[k].L(depth) * l;
	    mcell.nodes[j].N += projected_nodes[k].N(depth) * l;
	  }
	  
	}
	
	mesh.cells.push_back(mcell);
	refined_offset += refinedthickness;
      }
    }
  }

  void print(FILE *fp) const
  {
    fprintf(fp, "%f %d\n", thickness, (int)parameterset::NPARAMETERS);
    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      fprintf(fp, "%d ", (int)order[k]);
      
      for (size_t i = 0; i <= order[k]; i ++) {
	nodes[i].write_parameter(k, fp);
      }
      
      fprintf(fp, "\n");
    }
  }

  bool read(FILE *fp)
  {
    double fthickness;
    int fparameters;
    if (fscanf(fp, "%lf %d\n", &fthickness, &fparameters) != 2) {
      ERROR("Failed to reader header");
      return false;
    }

    if (fparameters != parameterset::NPARAMETERS) {
      ERROR("Invalid n parameters. (%d != %d)", fparameters, parameterset::NPARAMETERS);
      return false;
    }
    
    thickness = fthickness;

    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      int o;
      if (fscanf(fp, "%d", &o) != 1) {
	ERROR("Failed to read order");
	return false;
      }

      if (o < 0 || o > (int)maxorder) {
	ERROR("Order out of range");
	return false;
      }
      
      order[k] = o;
      for (size_t i = 0; i <= order[k]; i ++) {
	if (!nodes[i].read_parameter(k, fp)) {
	  ERROR("Failed to read node parameter");
	  return false;
	}
      }

    }
    
    return true;
  }

  bool save(FILE *fp) const
  {
    fprintf(fp, "%15.9f %d\n", (double)thickness, (int)parameterset::NPARAMETERS);

    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      fprintf(fp, "%d ", (int)order[k]);
      
      for (size_t i = 0; i <= order[k]; i ++) {
	
	if (!nodes[i].write_parameter(k, fp)) {
	  ERROR("Failed to save node parameter");
	  return false;
	}

      }

      fprintf(fp, "\n");
    }

    return true;
  }

  int encode_size()
  {
    int size = sizeof(int) + sizeof(double);

    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {
      size += sizeof(int);

      for (size_t i = 0; i <= order[k]; i ++) {
	size += nodes[i].encode_parameter_size(k);
      }
    }

    return size;
  }
  
  int encode(char *buffer, int &buffer_offset, int buffer_size)
  {
    int e;
    int s = 0;
    int n = (int)parameterset::NPARAMETERS;

    e = ::encode<int>(n, buffer, buffer_offset, buffer_size);
    if (e < 0) {
      ERROR("Failed to encode nparameters");
      return -1;
    }
    s += e;

    e = ::encode<double>((double)thickness, buffer, buffer_offset, buffer_size);
    if (e < 0) {
      ERROR("Failed to encode thickness");
      return -1;
    }
    s += e;

    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      int o = order[k];

      e = ::encode<int>(o, buffer, buffer_offset, buffer_size);
      if (e < 0) {
	ERROR("Failed to encode order");
	return -1;
      }
      s += e;
      
      for (size_t i = 0; i <= order[k]; i ++) {
	e = nodes[i].encode_parameter(k, buffer, buffer_offset, buffer_size);
	if (e < 0) {
	  ERROR("Failed to encode node %d", i);
	  return -1;
	}
	s += e;
      }
    }

    return s;
  }

  int decode(const char *buffer, int &buffer_offset, int buffer_size)
  {
    int n;

    if (::decode<int>(n, buffer, buffer_offset, buffer_size) < 0) {
      ERROR("Failed to decode no. parameters");
      return -1;
    }

    if (n != (int)parameterset::NPARAMETERS) {
      ERROR("Invalid no. parameters");
      return -1;
    }

    double t;
    if (::decode<double>(t, buffer, buffer_offset, buffer_size) < 0) {
      ERROR("Failed to encode thickness");
      return -1;
    }

    if (t <= 0.0) {
      ERROR("Invalid thickness");
      return -1;
    }

    thickness = t;

    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      int o;

      if (::decode<int>(o, buffer, buffer_offset, buffer_size) < 0) {
	ERROR("Failed to encode order");
	return -1;
      }

      if (o < 0 || o > (int)maxorder) {
	ERROR("Order out of range");
	return -1;
      }

      order[k] = o;
      
      for (size_t i = 0; i <= order[k]; i ++) {
	if (nodes[i].decode_parameter(k, buffer, buffer_offset, buffer_size) < 0) {
	  ERROR("Failed to encode node %d", i);
	  return -1;
	}
      }
    }
    return 0;
  }
    
  MeshParameter<real> lower_boundary(real offset) const
  {
    // Copy last node for each order to projected nodes
    
    for (size_t i = 0; i < parameterset::NPARAMETERS; i ++) {
      projected_nodes[maxorder] = nodes[order[i]];
    }

    MeshParameter<real> r;

    real depth = offset + thickness;
    
    r.rho = projected_nodes[maxorder].rho(depth);
    r.A = projected_nodes[maxorder].A(depth);
    r.C = projected_nodes[maxorder].C(depth);
    r.F = projected_nodes[maxorder].F(depth);
    r.L = projected_nodes[maxorder].L(depth);
    r.N = projected_nodes[maxorder].N(depth);
    
    return r;
  }

  parameterset interpolate(Mesh<real, maxorder> &mesh, real xi)
  {
    parameterset p;

    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      p[k] = 0.0;

      for (size_t i = 0; i <= order[k]; i ++) {

	p[k] += nodes[i][k] * mesh.quadrature[order[k]]->cardinal[i].value(xi);

      }
    }

    return p;
  }

  parameterset interpolate_gradient(Mesh<real, maxorder> &mesh,
				    real xi,
				    real depth,
				    size_t moffset,
				    Spec1DMatrix<real> &jacobian)
  {
    parameterset p;
    size_t poffset = 0;
    
    for (size_t k = 0; k < parameterset::NPARAMETERS; k ++) {

      p[k] = 0.0;

      for (size_t i = 0; i <= order[k]; i ++) {

	real w = mesh.quadrature[order[k]]->cardinal[i].value(xi);
	p[k] += nodes[i][k] * w;

	jacobian(poffset + i, 6*(moffset) + 0) += w * nodes[i].drho(k, depth);
	jacobian(poffset + i, 6*(moffset) + 1) += w * nodes[i].dA(k, depth);
	jacobian(poffset + i, 6*(moffset) + 2) += w * nodes[i].dC(k, depth);
	jacobian(poffset + i, 6*(moffset) + 3) += w * nodes[i].dF(k, depth);
	jacobian(poffset + i, 6*(moffset) + 4) += w * nodes[i].dL(k, depth);
	jacobian(poffset + i, 6*(moffset) + 5) += w * nodes[i].dN(k, depth);
      }

      poffset += (order[k] + 1);
    }

    return p;
  }
    
  real thickness;

  size_t order[parameterset::NPARAMETERS];
  // std::array<parameterset, maxorder + 1> nodes;
  // mutable std::array<parameterset, maxorder + 1> projected_nodes;
				    
  parameterset nodes[maxorder + 1];
  mutable parameterset projected_nodes[maxorder + 1];
  
};

#endif // cell_hpp
