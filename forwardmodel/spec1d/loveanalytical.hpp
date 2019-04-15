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
#ifndef loveanalytical_hpp
#define loveanalytical_hpp

template
<
  typename real
>
real loveanalytical(real omega, real c)
{
  static constexpr double H = 10.0e3;
  static constexpr double B1 = 3.0e3;
  static constexpr double B2 = 5.0e3;
  static constexpr double u1 = 2.8e3 * B1 * B1;
  static constexpr double u2 = 3.2e3 * B2 * B2;

  return

    tan(omega * H * sqrt(1.0/(B1 * B1) - 1.0/(c * c))) -

    (u2 * (sqrt(1.0/(c * c) - 1.0/(B2 * B2))))/
    (u1 * (sqrt(1.0/(B1 * B1) - 1.0/(c * c))));
  
}

#endif // loveanalytical_hpp

