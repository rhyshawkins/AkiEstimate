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
#ifndef parameterset_hpp
#define parameterset_hpp

#include <array>

#include "encodedecode.hpp"

template
<
  typename real,
  size_t set_size
>
class ParameterSet : public std::array<real, set_size> {
public:

  static constexpr size_t NPARAMETERS = set_size;

  ParameterSet()
  {
  }

  virtual ~ParameterSet()
  {
  }

  virtual real rho(real depth) const = 0;
  
  virtual real A(real depth) const = 0;
  
  virtual real C(real depth) const = 0;
  
  virtual real F(real depth) const = 0;
  
  virtual real L(real depth) const = 0;
  
  virtual real N(real depth) const = 0;

  virtual real drho(size_t i, real depth) const = 0;
  
  virtual real dA(size_t i, real depth) const = 0;
  
  virtual real dC(size_t i, real depth) const = 0;
  
  virtual real dF(size_t i, real depth) const = 0;
  
  virtual real dL(size_t i, real depth) const = 0;
  
  virtual real dN(size_t i, real depth) const = 0;
  
  virtual bool read(FILE *fp)
  {
    for (auto &r : *this) {
      double t;
      if (fscanf(fp, "%lf", &t) != 1) {
	return false;
      }
      r = t;      
    }

    return true;
  }

  virtual bool read_parameter(int k, FILE *fp)
  {
    if (k < 0 || k >= (int)set_size) {
      return false;;
    }

    double t;

    if (fscanf(fp, "%lf", &t) != 1) {
      return false;
    }

    (*this)[k] = t;
    return true;
  }

  virtual bool write(FILE *fp) const
  {
    for (auto &r : *this) {
      fprintf(fp, "%15.9f ", (double)r);
    }
    fprintf(fp, "\n");
    return true;
  }

  virtual bool write_parameter(int k, FILE *fp) const
  {
    if (k < 0 || k >= (int)set_size) {
      return false;;
    }

    fprintf(fp, "%15.9f", (double)((*this)[k]));
    return true;
  }

  virtual size_t encode_size() const
  {
    return sizeof(real) * set_size;
  }

  virtual size_t encode_parameter_size(int k) const
  {
    return sizeof(real);
  }

  int encode(char *buffer, int &buffer_offset, int buffer_size) const
  {
    for (auto &r : *this) {
      if (::encode<real>(r, buffer, buffer_offset, buffer_size) < 0) {
	return -1;
      }
    }

    return encode_size();
  }

  int encode_parameter(int k, char *buffer, int &buffer_offset, int buffer_size) const
  {
    if (k < 0 || k >= (int)set_size) {
      return -1;
    }
    
    if (::encode<real>((*this)[k], buffer, buffer_offset, buffer_size) < 0) {
      return -1;
    }

    return encode_parameter_size(k);
  }

  int decode(const char *buffer, int &buffer_offset, int buffer_size)
  {
    for (auto &r : *this) {
      if (::decode<real>(r, buffer, buffer_offset, buffer_size) < 0) {
	return -1;
      }
    }

    return encode_size();
  }
    
  int decode_parameter(int k, const char *buffer, int &buffer_offset, int buffer_size)
  {
    if (k < 0 || k >= (int)set_size) {
      return -1;
    }
    
    if (::decode<real>((*this)[k], buffer, buffer_offset, buffer_size) < 0) {
      return -1;
    }

    return encode_parameter_size(k);
  }
};

#endif // parameterset_hpp
