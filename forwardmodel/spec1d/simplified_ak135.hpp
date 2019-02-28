#pragma once
#ifndef simplified_ak135_hpp
#define simplified_ak135_hpp

#include <vector>

struct simplified_ak135entry {
  double depth;   // km
  double density; // Mg/km3
  double Vp;      // km/s
  double Vs;      // km/s
};

extern const std::vector<simplified_ak135entry> simplified_ak135_table;

double simplified_ak135_density(double depth);
double simplified_ak135_Vp(double depth);
double simplified_ak135_Vs(double depth);

template
<typename real>
struct simplified_ak135_1dreference {

  static real interpolate_rho(real depth)
  {
    return simplified_ak135_density(depth);
  }

  static real interpolate_Vp(real depth)
  {
    return simplified_ak135_Vp(depth);
  }
  
  static real interpolate_Vs(real depth)
  {
    return simplified_ak135_Vs(depth);
  }
};

#endif // simplified_ak135_hpp
