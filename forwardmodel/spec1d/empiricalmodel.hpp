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

