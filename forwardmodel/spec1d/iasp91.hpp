#pragma once
#ifndef iasp91_hpp
#define iasp91_hpp

#include <vector>

struct iasp91entry {
  double depth;   // km
  double density; // Mg/km3
  double Vp;      // km/s
  double Vs;      // km/s
};

extern const std::vector<iasp91entry> iasp91_table;

double iasp91_density(double depth);
double iasp91_Vp(double depth);
double iasp91_Vs(double depth);


#endif // iasp91_hpp
