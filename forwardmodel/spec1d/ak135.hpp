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
