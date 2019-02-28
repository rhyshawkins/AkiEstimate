#pragma once
#ifndef encodedecode_hpp
#define encodedecode_hpp

#include <string.h>

#include "logging.hpp"

template
<
  typename T
>
int encode(const T &t, char *buffer, int &offset, int buffer_size)
{
  if ((int)(offset + sizeof(T)) > buffer_size) {
    ERROR("Offset out of range: %d + %d > %d", offset, (int)sizeof(T), buffer_size);
    return -1;
  }

  memcpy(buffer + offset, &t, sizeof(T));
  offset += sizeof(T);
  return sizeof(T);
}

template
<
  typename T
>
int decode(T &t, const char *buffer, int &offset, int buffer_size)
{
  if ((int)(offset + sizeof(T)) > buffer_size) {
    ERROR("Offset out of range: %d + %d > %d", offset, (int)sizeof(T), buffer_size);
    return -1;
  }

  memcpy(&t, buffer + offset, sizeof(T));
  offset += sizeof(T);
  return sizeof(T);
}

#endif // encodedecode_hpp

