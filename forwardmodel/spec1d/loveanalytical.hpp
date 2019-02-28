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

