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
