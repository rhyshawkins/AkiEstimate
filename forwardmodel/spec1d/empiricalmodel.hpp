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
#ifndef empiricalmodel_hpp
#define empiricalmodel_hpp


template
<
  typename real
> class BrocherEmpiricalModel {
public:

  static void compute(const real &vs,
		      real &rho,
		      real &vp)
  {
    //
    // Empirical model from Brocher, 2005, BSSA
    // Modified for m/s
    //
    vp = 0.9409e3 + vs*(2.0947 + vs*(-0.8206e-3 + vs*(0.2683e-6 - 0.0251e-9*vs)));
    rho = vp*(1.6612 + vp*(-0.4721e-3 + vp*(0.0671e-6 + vp*(-0.0043e-9 + 0.000106e-12*vp))));


  }

  static void compute_gradient(const real &vs,
			       real &rho,
			       real &vp,
			       real &drhodvs,
			       real &dvpdvs)
  {
    //
    // Empirical model from Brocher, 2005, BSSA
    // Modified for m/s
    //
    vp = 0.9409e3 + vs*(2.0947 + vs*(-0.8206e-3 + vs*(0.2683e-6 - 0.0251e-9*vs)));
    dvpdvs = 2.0947 + vs*(-1.6412e-3 + vs*(0.8049e-6 - 0.1004e-9*vs));
			  
    rho = vp*(1.6612 + vp*(-0.4721e-3 + vp*(0.0671e-6 + vp*(-0.0043e-9 + 0.000106e-12*vp))));
    real drhodvp = 1.6612 + vp*(-0.9442e-3 + vp*(0.2013e-6 + vp*(-0.0172e-9 + 0.000530e-12*vp)));
    drhodvs = drhodvp * dvpdvs;
  }
};


#endif // empiricalmodel_hpp

