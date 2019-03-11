#pragma once
#ifndef dispersion_hpp
#define dispersion_hpp

#include <fftw3.h>

class DispersionData {
public:

  typedef enum {
    MAP_TYPE_CAUSAL,
    MAP_TYPE_ACAUSAL,
    MAP_TYPE_PRODUCT
  } te_map_t;
  
  DispersionData(double _fmin, double _fmax) :
    fmin(_fmin),
    fmax(_fmax),
    fullspec(NULL),
    fullenv(NULL),
    sigma_phase(0.1e03)
  {
  }

  bool load(const char *filename)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to open %s for reading\n", filename);
      return false;
    }

    if (fscanf(fp, "%lf %lf %lf %lf %lf\n",
               &lon1, &lat1, &lon2, &lat2, &distkm) != 5) {
      fprintf(stderr, "error: failed to parse line 1\n");
      return false;
    }

    if (fscanf(fp, "%lf %d %lf %lf %d\n",
	       &samplerate, &daycount, &asnr, &csnr, &samples) != 5) {
      fprintf(stderr, "error: failed to parse line 2\n");
      return false;
    }

    freq.resize(samples);
    sreal.resize(samples);
    simag.resize(samples);
    nreal.resize(samples);
    nimag.resize(samples);

    ffirst = samples;
    flast = 0;
    for (int i = 0; i < samples; i ++) {

      
      if (fscanf(fp, "%lf %lf %lf %lf %lf\n",
		 &freq[i],
		 &sreal[i], &simag[i],
		 &nreal[i], &nimag[i]) != 5) {
	fprintf(stderr, "error: failed to read spectrum\n");
	return false;
      }

      if (freq[i] >= fmin && i < ffirst) {
	ffirst = i;
      }

      if (freq[i] <= fmax && i > flast) {
	flast = i;
      }

    }

    fclose(fp);


    target_phase.resize(freq.size());
    target_error.resize(freq.size());

    predicted_phase.resize(freq.size());
    predicted_group.resize(freq.size());
    
    return true;
  }

  bool load_phase(const char *filename)
  {
    FILE *fp;

    fp = fopen(filename, "r");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to open phase file: %s\n", filename);
      return false;
    }

    while (!feof(fp)) {

      double lf, lc;
      int sign, offset;
      double le;
      
      if (fscanf(fp, "%lf %lf %d %d %lf\n", &lf, &lc, &sign, &offset, &le) != 5) {
	if (feof(fp)) {
	  break;
	} else {

	  fprintf(stderr, "error: failed to parse line\n");
	  return false;

	}
      }

      fpoints.push_back(lf);
      cpoints.push_back(lc * 1.0e3); // Optimization in metres
      epoints.push_back(le * 1.0e3);

    }

    fclose(fp);
    return true;
  }

  bool save_predictions(const char *filename)
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return false;
    }

    int N = freq.size();
    for (int i = 0; i < N; i ++) {
      fprintf(fp, "%15.9f %15.9f %15.9f %15.9f\n",
	      freq[i],
	      predicted_group[i],
	      predicted_phase[i],
	      target_phase[i]);
    }

    fclose(fp);
    return true;
  }


  double lerp(double f, double &e)
  {
    e = 0.0;
    
    if (f < fpoints[0]) {
      return 0.0;
    }

    if (f > fpoints[fpoints.size() - 1]) {
      return 0.0;
    }

    size_t i = 1;
    while (i < fpoints.size()) {

      if (f >= fpoints[i - 1] && f <= fpoints[i]) {

	double alpha = (f - fpoints[i - 1])/(fpoints[i] - fpoints[i - 1]);

	double c = cpoints[i - 1]*(1.0 - alpha) + cpoints[i]*alpha;
	e = epoints[i - 1]*(1.0 - alpha) + epoints[i]*alpha;
	
	return c;
      }

      i ++;
    }

    // Shouldn't get here
    return 0.0;
  }

  void initialise_target()
  {
    for (auto &c : target_phase) {
      c = 0.0;
    }
    for (auto &e : target_error) {
      e = 0.0;
    }

    bool started = false;
    
    for (int i = ffirst; i <= flast; i ++) {

      double e;
      double c = lerp(freq[i], e);
      if (c == 0.0) {
	if (started) {
	  flast = i - 1;
	  break;
	} else {
	  ffirst = i + 1;
	}
      } else {
	started = true;
	target_phase[i] = c;
	target_error[i] = e;
      }
    }

    //
    // Sanity
    //
    for (int i = ffirst; i <= flast; i ++) {
      if (target_phase[i] == 0.0) {
	fprintf(stderr, "target phase not set\n");
	throw std::exception();
      }
    }
      
  }

  
  double lon1, lat1;
  double lon2, lat2;
  double distkm;

  double samplerate;
  int daycount;
  double asnr, csnr;
  int samples;

  std::vector<double> freq, sreal, simag, nreal, nimag;
  double fmin;
  double fmax;

  int ffirst, flast;
  
  fftw_plan plan;
  fftw_complex *fullspec;
  fftw_complex *fullenv;

  std::vector<double> predicted_group;
  std::vector<double> predicted_phase;

  double sigma_phase;
  std::vector<double> target_phase;
  std::vector<double> target_error;

  std::vector<double> fpoints, cpoints, epoints;
};

#endif // dispersion_hpp
    

  
