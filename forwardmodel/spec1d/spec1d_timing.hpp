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
#ifndef spec1d_timing_hpp
#define spec1d_timing_hpp

#include <time.h>

class Spec1DTimer {
public:

  Spec1DTimer()
  {
  }

  void start()
  {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
  }

  double stop()
  {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);

    return (double)(stop_time.tv_sec - start_time.tv_sec) * 1.0e6 +
      (double)(stop_time.tv_nsec - start_time.tv_nsec)/1.0e3;
  }

private:

  timespec start_time, stop_time;
  
};

#endif // spec1d_timing_hpp
  
  
