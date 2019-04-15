//
//    AkiEstimate : A method for the joint estimation of Love and Rayleigh surface wave
//    dispersion from ambient noise cross-correlations.
//
//      Hawkins R. and Sambridge M., "An adjoint technique for estimation of interstation phase
//    and group dispersion from ambient noise cross-correlations", BSSA, 2019
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
#ifndef dispersion_hpp
#define dispersion_hpp

#include <fftw3.h>
#include <gsl/gsl_sf_bessel.h>

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
    env_signal(NULL),
    env_spectrum(NULL),
    noise_sigma(1.0)
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
    ncfreal.resize(samples);
    ncfimag.resize(samples);

    predicted_k.resize(samples);
    predicted_group.resize(samples);
    predicted_phase.resize(samples);
    predicted_bessel.resize(samples);
    predicted_envelope.resize(samples);
    predicted_realspec.resize(samples);

    ffirst = samples;
    flast = 0;
    for (int i = 0; i < samples; i ++) {

      
      if (fscanf(fp,
		 "%lf %lf %lf %lf %lf\n",
		 &freq[i],
		 &sreal[i], &simag[i],
		 &ncfreal[i], &ncfimag[i]) != 5) {
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
      fprintf(fp, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n",
	      freq[i],
	      predicted_group[i],
	      predicted_phase[i],
	      predicted_bessel[i],
	      predicted_envelope[i],
	      predicted_realspec[i]);
    }

    fclose(fp);
    return true;
  }

  bool compute_ftan_envelope(double cfreq, double sigma, double *acausal, double *causal)
  {
    //
    // Effectively computes | H ( F^-1( G(cfreq, sigma) * (sreal + j*simag)) ) | where H is the
    // Hilbert transform, F^-1 the inverse Fourier transform, G a Gaussian filter.
    //
    int N2 = samples - 1;
    int N = N2 * 2;
    if (fullspec == NULL) {

      fullspec = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
      fullenv = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

      plan = fftw_plan_dft_1d(N, fullspec, fullenv, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    //
    // Set -ve frequencies to zero and positive to twice the filtered signal
    //
    for (int i = 0; i < samples; i ++) {

      if (i > 0) {
	fullspec[i][0] = 0.0;
	fullspec[i][1] = 0.0;
      }

      double df = freq[i] - cfreq;
      double g = 2.0*exp(-(df*df)/(2.0 * sigma*sigma));
      fullspec[i + samples][0] = g * sreal[i];
      fullspec[i + samples][1] = g * simag[i];
    }

    //
    // Inverse FFT
    //
    fftw_execute(plan);

    //
    // Absolute value gives envelope
    //
    for (int i = 0; i < N2; i ++) {
      acausal[i] = env(fullenv[N - 1 - i][0], fullenv[N - 1 - i][1]);
      causal[i] = env(fullenv[i][0], fullenv[i][1]);
    }

    return true;
  }

  double env(double r, double i)
  {
    return sqrt(r*r + i*i);
  }

  void estimate_sigma(double fthresh)
  {
    // int meann = 0;
    // double mean = 0.0;
    // double var = 0.0;
    
    // for (int i = 0; i < (int)ncfreal.size(); i ++) {
    //   if (freq[i] >= fthresh) {
    // 	meann ++;

    // 	double delta = ncfreal[i] - mean;
    // 	mean += delta/(double)meann;
    // 	var += delta*(ncfreal[i] - mean);
    //   }
    // }

    double maxA = 0.0;
    for (int i = ffirst; i <= flast; i ++) {
      double A = fabs(ncfreal[i]);
      if (A > maxA) {
	maxA = A;
      }
    }
      
    noise_sigma = maxA * 0.05;
  }
  
  bool build_time_energy_map(double sigma,
			     te_map_t map_type,
			     double vmin,
			     double vmax,
			     double max_deltav)
  {
    int N2 = samples - 1;
    int N = N2 * 2;

    //
    // Create
    //
    time_energy_map.resize(freq.size());
    amplitude_max.resize(freq.size());
    time_max.resize(freq.size());
    time.resize(N);

    predicted_group.resize(freq.size());
    predicted_phase.resize(freq.size());
    predicted_bessel.resize(freq.size());
    predicted_envelope.resize(freq.size());
    predicted_realspec.resize(freq.size());
    
    //
    // Initialize
    //
    for (auto &t : time_energy_map) {
      t = nullptr;
    }
    for (auto &t : time_max) {
      t = -1.0;
    }

    double tmin = distkm/vmax;
    double tmax = distkm/vmin;

    int itmin = N;
    int itmax = 0;
    
    for (int j = 0; j < N; j ++) {
      time[j] = ((double)j + 0.5)/samplerate;

      if (time[j] >= tmin && j < itmin) {
	itmin = j;
      }
      if (time[j] <= tmax && j > itmax) {
	itmax = j;
      }
    }

    double *causal = new double[N];
    double *acausal = new double[N];

    double df = freq[1] - freq[0];
    double dt = 1.0/samplerate;

    int last_maxj = -1;

    max_amplitude = 0.0;

    for (int i = ffirst; i <= flast; i ++) {

      time_energy_map[i] = new double[N];

      if (!compute_ftan_envelope(freq[i], sigma, acausal, causal)) {
	return false;
      }

      switch (map_type) {
      case MAP_TYPE_CAUSAL:
	for (int j = 0; j < N; j ++) {
	  time_energy_map[i][j] = causal[j];
	}
	break;

      case MAP_TYPE_ACAUSAL:
	for (int j = 0; j < N; j ++) {
	  time_energy_map[i][j] = acausal[j];
	}
	break;

      case MAP_TYPE_PRODUCT:
	for (int j = 0; j < N; j ++) {
	  time_energy_map[i][j] = causal[j] * acausal[j];
	}
	break;
	
      default:
	fprintf(stderr, "error: invalid map type\n");
	return false;
      }

      int maxj = -1;
      if (last_maxj > 0 && max_deltav > 0.0) {

	double deltaU = max_deltav * df;
	double lastt = time[last_maxj];
	double lastU = distkm/lastt;
	double deltat = lastt - distkm/(lastU + deltaU);

	int mvitmin = last_maxj - (int)ceil(deltat/dt);
	if (mvitmin < itmin) {
	  mvitmin = itmin;
	}

	double maxe = 0.0;
	for (int j = mvitmin; j <= itmax; j ++) {
	  if (time_energy_map[i][j] > maxe) {
	    maxe = time_energy_map[i][j];
	    maxj = j;
	  }
	}

      } else {
	double maxe = 0.0;
	for (int j = itmin; j <= itmax; j ++) {
	  if (time_energy_map[i][j] > maxe) {
	    maxe = time_energy_map[i][j];
	    maxj = j;
	  }
	}
      }

      last_maxj = maxj;
	
      time_max[i] = time[maxj];
      amplitude_max[i] = time_energy_map[i][maxj];
      if (amplitude_max[i] > max_amplitude) {
	max_amplitude = amplitude_max[i];
      }
    }

    delete causal;
    delete acausal;

    return true;
  }

  bool save_max_time(const char *filename)
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to create %s\n", filename);
      return false;
    }

    int N = freq.size();
    for (int i = 0; i < N; i ++) {
      double t = time_max[i];
      double U = 0.0;
      if (t > 0.0) {
	U = distkm/t;
      } else {
	t = 0.0;
      }
      fprintf(fp, "%15.9f %16.9f %16.9f %16.9e\n", freq[i], t, U * 1.0e3, amplitude_max[i]);
    }

    fclose(fp);
    return true;
  }

  bool save_time_energy(const char *filename)
  {
    int N2 = samples - 1;
    int N = N2 * 2;

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to create %s\n", filename);
      return false;
    }

    for (auto &p : time_energy_map) {

      if (p == nullptr) {
	//
	// Output zero row
	//
	for (int i = 0; i < N; i ++) {
	  fprintf(fp, "0.0 ");
	}
	
      } else {

	for (int i = 0; i < N; i ++) {
	  fprintf(fp, "%16.9e ", p[i]);
	}
      }

      fprintf(fp, "\n");
    }

    fclose(fp);

    return true;
  }

  void compute_bessel()
  {
    for (int i = ffirst; i <= flast; i ++) {

      if (predicted_phase[i] > 0.0) {

	double k = 2.0 * M_PI * freq[i]/predicted_phase[i];
	predicted_bessel[i] = gsl_sf_bessel_J0(k * distkm * 1.0e3);

      }

    }
  }

  void set_unit_envelope()
  {
    for (auto &e : predicted_envelope) {
      e = 1.0;
    }
  }
  
  void compute_envelope(double gaussian_smooth_sigma)
  {
    //
    // The envelope is computed taking the envelope of the real part of
    // the spectrum.  
    //


    //
    // Compute the envelope of the real part of the spectrum. Once off
    // at start. Little trick here, spectrum is saved as power 2 + 1, eg 4097
    // samples, so we chop one off to get a power of 2 fft.
    //
    if (env_signal == NULL) {
      
      env_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (samples - 1));
      env_spectrum = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (samples - 1));

      env_fplan =
	fftw_plan_dft_1d(samples - 1, env_signal, env_spectrum, FFTW_FORWARD, FFTW_ESTIMATE);
      env_bplan =
	fftw_plan_dft_1d(samples - 1, env_spectrum, env_signal, FFTW_BACKWARD, FFTW_ESTIMATE);
      
    }

    // Copy signal
    for (int i = 1; i < samples; i ++) {
      env_signal[i - 1][0] = ncfreal[i];
      env_signal[i - 1][1] = 0.0;
    }

    // Execute
    fftw_execute(env_fplan);
    
    // Zero negative frequencies and double positive
    int half = (samples - 1)/2;
    for (int i = 0; i < half; i ++) {
      env_spectrum[i][0] *= 2.0;
      env_spectrum[i][1] *= 2.0;
      env_spectrum[half + i][0] = 0.0;
      env_spectrum[half + i][1] = 0.0;
    }
      
    // Inverse transform
    fftw_execute(env_bplan);

    // Absolute value is envelope
    for (int i = 1; i < samples; i ++) {
      predicted_envelope[i] = env(env_signal[i - 1][0]/(double)(samples - 1), env_signal[i - 1][1]/(double)(samples - 1));
    }
    predicted_envelope[0] = 0.0;

    // Smoothing filter
    if (gaussian_smooth_sigma > 0.0) {
      int hn = (int)ceil(3.0 * gaussian_smooth_sigma);
      if (hn == 0) {
	hn ++;
      }
      
      // Build filter
      int fn = 2*hn - 1;
      double *filter = new double[fn];
      for (int i = 0; i < fn; i ++) {
	filter[i] = 0.0;
      }
      
      double sfilter = 0.0;
      for (int i = 0; i < hn; i ++) {
	double dx = (double)i;
	double g = exp(-(dx * dx)/(2.0 * gaussian_smooth_sigma * gaussian_smooth_sigma));
	
	if (i == 0) {
	  filter[hn - 1] = g;
	  sfilter += g;
	} else {
	  filter[hn - 1 + i] = g;
	  filter[hn - 1 - i] = g;
	  sfilter += 2.0*g;
	}
      }
      
      // Normalise
      for (int i = 0; i < fn; i ++) {
	filter[i] /= sfilter;
      }
      
      // Do filter
      for (int i = 0; i < samples; i ++) {
	
	double s = 0.0;
	for (int j = 0; j < fn; j ++) {
	  
	  int si = i + j - hn;
	  if (si < 0) {
	    si = 0;
	  } else if (si >= samples) {
	    si = samples - 1;
	  }
	  
	  s += predicted_envelope[si] * filter[j];
	}
	
	// Use real spec array as a temporary
	predicted_realspec[i] = s;
      }
      
      // Copy back
      for (int i = 0; i < samples; i ++) {
	predicted_envelope[i] = predicted_realspec[i];
      }
      
      delete [] filter;
	
    }

  }
  
  double lon1, lat1;
  double lon2, lat2;
  double distkm;

  double samplerate;
  int daycount;
  double asnr, csnr;
  int samples;

  std::vector<double> freq, sreal, simag, ncfreal, ncfimag;
  double fmin;
  double fmax;

  int ffirst, flast;
  
  fftw_plan plan;
  fftw_complex *fullspec;
  fftw_complex *fullenv;

  fftw_plan env_fplan, env_bplan;
  fftw_complex *env_signal;
  fftw_complex *env_spectrum;

  std::vector<double *> time_energy_map;
  std::vector<double> amplitude_max;
  std::vector<double> time_max;
  std::vector<double> time;
  double max_amplitude;

  std::vector<double> predicted_k;
  std::vector<double> predicted_group;
  std::vector<double> predicted_phase;
  std::vector<double> predicted_bessel;
  std::vector<double> predicted_envelope;
  std::vector<double> predicted_realspec;

  double noise_sigma;
};

#endif // dispersion_hpp
    

  
