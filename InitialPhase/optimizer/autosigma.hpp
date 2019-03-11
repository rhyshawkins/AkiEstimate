#ifndef autosigma_hpp
#define autosigma_hpp

#include <math.h>

double autosigma(double distkm)
{
  constexpr double SLOPE = 0.441895439828;
  constexpr double INTERCEPT = -0.00729666600638;

  return SLOPE/sqrt(distkm) + INTERCEPT;
}

#endif // autosigma_hpp
