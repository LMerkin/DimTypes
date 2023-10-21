// vim:ts=2:et
//===========================================================================//
//                            "Tests/DimTests.cpp":                          //
//===========================================================================//
#include "DimTypes/DimTypes.hpp"
#include <cstdio>

namespace
{
  //=========================================================================//
  // Type Declarations:                                                      //
  //=========================================================================//
  DECLARE_DIMS(
    double,
    (Len,  m,   (km,  1000.0),  (AU, 1.495978706996262e+11)),
    (Time, sec, (day, 86400.0)),
    (Mass, kg)
  )
}

int main()
{
  using namespace std;
  printf("MaxHeight=%d\n\n", DimTypes::Bits::MaxHeight);

  //=========================================================================//
  // Astronomical Constants (from DE423):                                    //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // Speed of light:                                                         //
  //-------------------------------------------------------------------------//
  constexpr auto c = 299792.458_km / 1.0_sec;

  //-------------------------------------------------------------------------//
  // Heliocentric Gravitational Constant:                                    //
  //-------------------------------------------------------------------------//
  constexpr auto GMS  =
    2.959122082855911e-4 * IPow<3>(1.0_AU) / IPow<2>(1.0_day);

  //-------------------------------------------------------------------------//
  // Earth / Moon Mass Ratio: dimension-less:                                //
  //-------------------------------------------------------------------------//
  constexpr double EMRat = 81.30056941599857;

  //-------------------------------------------------------------------------//
  // Geocentric Gravitational Constant:                                      //
  //-------------------------------------------------------------------------//
  constexpr auto GME  =
    8.997011408268049e-10  / (1.0 + 1.0 / EMRat) *
    IPow<3>(1.0_AU) / IPow<2>(1.0_day);

  auto tonne   = Mass(1000.0);  // Assuming "kg" which is the FundUnit
  auto tn1     = CbRt(tonne);
  printf("CbRt(tonne) = %s\n", ToStr(tn1).data());

  auto tn2     = - tn1;
  auto tn3     = Abs(tn2);
  printf("tn3  = %s\n",   ToStr(tn3).data());

  auto AU1     = Len_AU(1.0);
  printf("AU1  = %s\n",   ToStr(AU1).data());
  auto AU2     = 1.0_AU;
  printf("AU2  = %s\n",   ToStr(AU2).data());
  auto AU3     = 1.0  * AU2;
  printf("AU3  = %s\n",   ToStr(AU3).data());
  auto AU4     = To_Len_km(AU1);
  printf("AU4  = %s\n",   ToStr(AU4).data());
  auto AU5     = To_Len_km(AU2);
  printf("AU5  = %s\n",   ToStr(AU5).data());

  printf("c    = %s\n",   ToStr(c)  .data());
  printf("GMS  = %s\n",   ToStr(GMS).data());
  constexpr auto kS     = SqRt (GMS);
  printf("kS   = %s\n",   ToStr(kS) .data());

  auto GMS1    = To_Time_sec(To_Len_km(GMS));
  printf("GMS1 = %s\n",   ToStr(GMS1).data());

  printf("GME  = %s\n",   ToStr(GME).data());
  constexpr auto GME1   = To_Time_sec(To_Len_km(GME));
  printf("GME1 = %s\n",   ToStr(GME1).data());
  constexpr auto kE     = SqRt (GME1);
  printf("kE   = %s\n",   ToStr(kE).data());

  auto x       = 10.0_km / 1.0_sec;
  auto y       = To_Len_AU(To_Time_day(x));
  auto z       = 1.0 / y;
  auto dl      = y * z;
  auto cmx     = c - x;

  printf("x    = %s\n",   ToStr(x)  .data());
  printf("     = %s\n",   ToStr(y)  .data());
  printf("1/x  = %s\n",   ToStr(z)  .data());
  printf("x/x  = %s\n",   ToStr(dl) .data());
  printf("c-x  = %s\n",   ToStr(cmx).data());
  return 0;
}
