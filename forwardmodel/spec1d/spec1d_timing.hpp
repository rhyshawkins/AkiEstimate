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
  
  
