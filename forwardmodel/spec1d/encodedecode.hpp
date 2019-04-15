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

