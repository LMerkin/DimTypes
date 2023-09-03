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
  /*
  DECLARE_DIMS(     Len,  Time,   Mass)

  DECLARE_DIM_UNITS(Len,  double, km,  AU)
  DECLARE_UNIT(     Len,  double, AU,  1.495978706996262e+8)

  DECLARE_DIM_UNITS(Time, double, sec, day)
  DECLARE_UNIT(     Time, double, day, 86400.0)

  DECLARE_DIM_UNITS(Mass, double, kg)

  DECLARE_DIM_STR
  */
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
  constexpr auto c = 299792.458 * Len_km / Time_sec;

  //-------------------------------------------------------------------------//
  // Heliocentric Gravitational Constant:                                    //
  //-------------------------------------------------------------------------//
  constexpr auto GMS  =
    2.959122082855911e-4 * IPow<3>(Len_AU) / IPow<2>(Time_day);

  //-------------------------------------------------------------------------//
  // Earth / Moon Mass Ratio: dimension-less:                                //
  //-------------------------------------------------------------------------//
  constexpr double EMRat = 81.30056941599857;

  //-------------------------------------------------------------------------//
  // Geocentric Gravitational Constant:                                      //
  //-------------------------------------------------------------------------//
  constexpr auto GME  =
    8.997011408268049e-10  / (1.0 + 1.0 / EMRat) *
    IPow<3>(Len_AU) / IPow<2>(Time_day);

  auto tonne = Mass_kg_T(1000.0);
  auto tn1   = CbRt(tonne);
  printf("CbRt(tonne) = %s\n", ToStr(tn1));

  auto tn2 = - tn1;
  auto tn3 = Abs(tn2);
  printf("tn3  = %s\n",   ToStr(tn3));

  auto AU1     = Len_AU_T(1.0);
  auto AU2     = 1.0 * Len_AU;
  auto AU3     = To_Len_km(AU1);
  auto AU4     = To_Len_km(AU2);
  printf("AU1  = %s\n",   ToStr(AU1));
  printf("AU2  = %s\n",   ToStr(AU2));
  printf("AU3  = %s\n",   ToStr(AU3));
  printf("AU4  = %s\n",   ToStr(AU4));

  printf("c    = %s\n",   ToStr(c));
  printf("GMS  = %s\n",   ToStr(GMS));
  printf("GME  = %s\n",   ToStr(GME));

  auto GMS1 = To_Time_sec(To_Len_km(GMS));
  auto GME1 = To_Time_sec(To_Len_km(GME));

  printf("GMS1 = %s\n",   ToStr(GMS1));
  printf("GME1 = %s\n",   ToStr(GME1));

  auto x       = 10.0 * Len_km / Time_sec;
  auto y       = To_Len_AU(To_Time_day(x));
  auto z       = 1.0 / y;
  auto dl      = y * z;
  auto cmx     = c - x;

  printf("x    = %s\n",   ToStr(x));
  printf("     = %s\n",   ToStr(y));
  printf("1/x  = %s\n",   ToStr(z));
  printf("x/x  = %s\n",   ToStr(dl));
  printf("c-x  = %s\n",   ToStr(cmx));
  return 0;
}
