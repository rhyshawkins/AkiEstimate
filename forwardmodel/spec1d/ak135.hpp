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
#ifndef ak135_hpp
#define ak135_hpp

#include <vector>

struct ak135entry {
  double depth;   // km
  double density; // Mg/km3
  double Vp;      // km/s
  double Vs;      // km/s
};

extern const std::vector<ak135entry> ak135_table;

double ak135_density(double depth);
double ak135_Vp(double depth);
double ak135_Vs(double depth);

template
<typename real>
struct ak135_1dreference {

  static real interpolate_rho(real depth)
  {
    return ak135_density(depth);
  }

  static real interpolate_Vp(real depth)
  {
    return ak135_Vp(depth);
  }
  
  static real interpolate_Vs(real depth)
  {
    return ak135_Vs(depth);
  }
};

#endif // ak135_hpp
