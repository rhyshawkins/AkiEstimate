#pragma once
#ifndef modelhistory_hpp
#define modelhistory_hpp

#include "model.hpp"
#include "encodedecode.hpp"
template
<
  typename real,
  typename parameterset,
  typename boundarycondition,
  size_t maxorder
>
class ModelHistoryWriter {
public:

  typedef Model<real, parameterset, boundarycondition, maxorder> model_t;
  
  ModelHistoryWriter(const char *filename, size_t _buffer_size = 1048576) :
    fp(fopen(filename, "w")),
    buffer(new char[_buffer_size]),
    buffer_size(_buffer_size),
    buffer_offset(0)
  {

    if (fp == NULL) {
      FATAL("Failed to create file.");
    }
    
  }

  ~ModelHistoryWriter()
  {
    flush_buffer();
    fclose(fp);

    delete [] buffer;
  }

  void add(model_t &model,
	   double likelihood,
	   int nhyper,
	   const real *hyper,
	   int npred,
	   const real *pred)
  {
    int size = model.encode_size();

    size += 2*sizeof(int) + (nhyper + npred + 1) * sizeof(double);

    if ((int)(size + sizeof(int)) > buffer_size) {
      FATAL("Buffer size too small");
    }

    if ((int)(buffer_offset + size + sizeof(int)) > buffer_size) {
      flush_buffer();
    }

    if (encode<int>(size, buffer, buffer_offset, buffer_size) < 0) {
      FATAL("Failed to encode size");
    }

    if (encode<double>(likelihood, buffer, buffer_offset, buffer_size) < 0) {
      FATAL("Failed to encode likelihood");
    }

    if (encode<int>(nhyper, buffer, buffer_offset, buffer_size) < 0) {
      FATAL("Failed to encode no. hyper parameters");
    }

    for (int i = 0; i < nhyper; i ++) {
      if (encode<double>(hyper[i], buffer, buffer_offset, buffer_size) < 0) {
	FATAL("Failed to encode hyper parameter");
      }
    }

    if (encode<int>(npred, buffer, buffer_offset, buffer_size) < 0) {
      FATAL("Failed to encode no. hyper parameters");
    }

    for (int i = 0; i < npred; i ++) {
      if (encode<double>(pred[i], buffer, buffer_offset, buffer_size) < 0) {
	FATAL("Failed to encode prediction");
      }
    }

    if (model.encode(buffer, buffer_offset, buffer_size) < 0) {
      FATAL("Failed to encode model");
    }
  }

private:

  void flush_buffer()
  {
    if (buffer_offset > 0) {

      if (fwrite(buffer, buffer_offset, 1, fp) != 1) {
	FATAL("Failed to write to file");
      }

      buffer_offset = 0;
    }
  }
  
  FILE *fp;
  char *buffer;
  int buffer_size;
  int buffer_offset;
  
};

template
<
  typename real,
  typename parameterset,
  typename boundarycondition,
  size_t maxorder
>
class ModelHistoryReader {
public:

  typedef Model<real, parameterset, boundarycondition, maxorder> model_t;

  ModelHistoryReader(const char *filename, int _buffer_size = 65536) :
    fp(fopen(filename, "r")),
    buffer(new char[_buffer_size]),
    buffer_size(_buffer_size)
  {
    if (fp == NULL) {
      FATAL("Failed to open file");
    }
  }

  ~ModelHistoryReader()
  {
    fclose(fp);
    delete [] buffer;
  }

  int next(model_t &model,
	   double &likelihood,
	   int maxhyper,
	   int &nhyper,
	   double *hyper,
	   int maxpred,
	   int &npred,
	   double *pred)
  {
    int size;

    if (fread(&size, sizeof(int), 1, fp) != 1) {
      if (feof(fp)) {
	return 0;
      } else {
	ERROR("Failed to read size of next model step");
	return -1;
      }
    }

    check_buffer_size(size);
    if (fread(buffer, size, 1, fp) != 1) {
      ERROR("Failed to read next %d bytes for next model", size);
      return -1;
    }

    int buffer_offset = 0;

    if (::decode<double>(likelihood, buffer, buffer_offset, size) < 0) {
      ERROR("Failed to decode likelihood");
      return -1;
    }

    if (::decode<int>(nhyper, buffer, buffer_offset, size) < 0) {
      ERROR("Failed to decode no. hyper parameters");
      return -1;
    }

    if (nhyper > maxhyper) {
      ERROR("Too many hyper parameters");
      return -1;
    }

    for (int i = 0; i < nhyper; i ++) {
      if (::decode<double>(hyper[i], buffer, buffer_offset, size) < 0) {
	return -1;
      }
    }
    
    if (::decode<int>(npred, buffer, buffer_offset, size) < 0) {
      ERROR("Failed to decode no. pred parameters");
      return -1;
    }

    if (npred > maxpred) {
      ERROR("Too many pred parameters");
      return -1;
    }

    for (int i = 0; i < npred; i ++) {
      if (::decode<double>(pred[i], buffer, buffer_offset, size) < 0) {
	return -1;
      }
    }

    if (model.decode(buffer, buffer_offset, size) < 0) {
      ERROR("Failed to decode model");
      return -1;
    }

    return 1;
  }
	   

private:

  void check_buffer_size(int size)
  {
    int new_buffer_size = buffer_size;

    while (new_buffer_size < size) {
      new_buffer_size *= 2;
    }

    if (new_buffer_size != buffer_size) {
      delete [] buffer;
      buffer = new char[new_buffer_size];
      buffer_size = new_buffer_size;
    }
  }

  FILE *fp;
  char *buffer;
  int buffer_size;

};

#endif // modelhistory_hpp
