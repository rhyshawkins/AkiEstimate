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


#ifndef halfspace_hpp
#define halfspace_hpp

#include <stdio.h>
#include <string>

template
<
  typename real,
  size_t np
>
class Halfspace {
public:

  Halfspace()
  {
  }

  virtual ~Halfspace()
  {
  }

  bool read(FILE *fp)
  {
    int nparameters;
    double thickness;
    
    if (fscanf(fp, "%d %lf\n", &nparameters, &thickness) != 2) {
      ERROR("Failed to read cell parameters");
      return false;
    }

    if (nparameters != (int)np || thickness != 0.0) {
      ERROR("Expected terminating cell parameters, got %d %f", nparameters, thickness);
      return false;
    }

    for (size_t i = 0; i < np; i ++) {
      double p;
      if (fscanf(fp, "%lf\n", &p) != 1) {
	ERROR("Failed to read parameters");
	return false;
      }

      parameters[i] = p;
    }

    return true;
  }
  
  bool save(FILE *fp) const
  {
    fprintf(fp, "%d 0.0\n", (int)np);
    for (size_t i = 0; i < np; i ++) {
      fprintf(fp, "%15.9f\n", (double)parameters[i]);
    }
    
    return true;
  }

  void print(FILE *fp) const
  {
    for (size_t i = 0; i < np; i ++) {
      fprintf(fp, "  %15.9f\n", (double)parameters[i]);
    }
  }

  int encode_size()
  {
    return np * sizeof(real);
  }
  
  int encode(char *buffer, int &buffer_offset, int buffer_size)
  {
    for (size_t i = 0; i < np; i ++) {
      if (::encode<real>(parameters[i], buffer, buffer_offset, buffer_size) < 0) {
	return -1;
      }
    }

    return np*sizeof(real);
  }

  int decode(const char *buffer, int &buffer_offset, int buffer_size)
  {
    for (size_t i = 0; i < np; i ++) {
      if (::decode<real>(parameters[0], buffer, buffer_offset, buffer_size) < 0) {
	return -1;
      }
    }
      
    return np*sizeof(real);
  }

  const real &operator[](int i) const
  {
    return parameters[i];
  }

  real &operator[](int i) 
  {
    return parameters[i];
  }

  static size_t nparameters()
  {
    return np;
  }
  
  static bool is_fixed()
  {
    return np == 0;
  }

  std::array<real, np> parameters;
  
};

#endif // halfspace_hpp
