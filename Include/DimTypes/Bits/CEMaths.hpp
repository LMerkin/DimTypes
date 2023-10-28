// vim:ts=2:et
//===========================================================================//
//                                 "CEMaths.hpp":                            //
//      Standard Mathematical Functions as "constexpr"s, for C++ <= 23       //
//===========================================================================//
#pragma  once
#include <cmath>
#include <cassert>
#include <cfloat>
#include <climits>
#include <complex>
#include <utility>
#include <cstdio>

namespace DimTypes::Bits::CEMaths
{
  //=========================================================================//
  // "IsComplex":                                                            //
  //=========================================================================//
  template<typename T>
  constexpr bool IsComplex = false;

  template<typename T>
  constexpr bool IsComplex<std::complex<T>> = true;

  //=========================================================================//
  // Consts:                                                                 //
  //=========================================================================//
  // Infinity, NaN:
  template<typename F> constexpr F Inf    = F(INFINITY);      // +oo
  template<typename F> constexpr F NaN    = F(NAN);           // NaN

  // Eps:
  template<typename F> F Eps;
  template<> constexpr float       Eps<float>       = FLT_EPSILON;
  template<> constexpr double      Eps<double>      = DBL_EPSILON;
  template<> constexpr long double Eps<long double> = LDBL_EPSILON;

  // Pi:
  template<typename F> constexpr F Pi     = F(M_PIl);         // Pi
  template<typename F> constexpr F TwoPi  = F(2.0L * M_PIl);  // 2*Pi
  template<typename F> constexpr F Pi_2   = F(M_PI_2l);       // Pi/2
  template<typename F> constexpr F Pi_4   = F(M_PI_4l);       // Pi/4
  template<typename F> constexpr F Pi_8   = Pi_4<F> / F(2.0); // Pi/8
  template<typename F> constexpr F Pi3_8  = Pi_4<F> * F(1.5); // 3/8 Pi

  // Ln2, Log2E
  template<typename F> constexpr F Ln2    = F(M_LN2l);        // ln(2)
  template<typename F> constexpr F Log2E  = F(M_LOG2El);      // log2(e)=1/Ln2

  // Cos and Sin of Pi/8:
  template<typename F> F CosPi8;
  template<> constexpr long double CosPi8<long double> =
    0.923879532511286756128L;
  template<> constexpr double      CosPi8<double>      =
    double(CosPi8<long double>);
  template<> constexpr float       CosPi8<float>       =
    float (CosPi8<long double>);

  template<typename F> F SinPi8;
  template<> constexpr long double SinPi8<long double> =
    0.382683432365089771728L;
  template<> constexpr double      SinPi8<double>      =
    double(SinPi8<long double>);
  template<> constexpr float       SinPi8<float>       =
    float (SinPi8<long double>);

  //=========================================================================//
  // Pade Approximants for Exp, Sin and Cos:                                 //
  //=========================================================================//
  // Different orders of Pade approximants are implemented, depending on the
  // actual floating-point type ("float", "double" or "long double"):
  //=========================================================================//
  // "ExpPade":                                                              //
  //=========================================================================//
  // The arg is assumed to be in [-0.5 .. 0.5]:
  //
  template<typename F>
  F ExpPade(F a_x);

  //-------------------------------------------------------------------------//
  // ExpPade<float>: 4th-Order Pade Approximant:                             //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-23) =~ 1.2e-7,
  // the abs precision of the 4th-order approximant is ~1e-10 (however, the
  // 3rd-order one would be slightly insuffient):
  //
  template<>
  constexpr float ExpPade<float>(float a_x)
  {
    assert(std::fabs(a_x) < 0.5f);

    constexpr float a1 = 0.5f;
    constexpr float a2 = float(3.0 /   28.0);
    constexpr float a3 = float(1.0 /   84.0);
    constexpr float a4 = float(1.0 / 1680.0);

    float  x2  = a_x * a_x;
    float  x3  = x2  * a_x;
    float  x4  = x2  * x2;

    float  t1  = a1  * a_x;
    float  t2  = a2  * x2;
    float  t3  = a3  * x3;
    float  t4  = a4  * x4;

    float  se  = t4  + t2  + 1.0f;
    float  so  = t3  + t1;
    float  res = (se + so) / (se - so);
    assert(res > 0.0f);
    return res;
  }

  //-------------------------------------------------------------------------//
  // ExpPade<double>: 6th-Order Pade Approximant:                            //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-52) =~ 2.2e-16,
  // the abs precision of the 6th-order approximation is ~4e-17:
  //
  template<>
  constexpr double ExpPade<double>(double a_x)
  {
    assert(std::fabs(a_x) < 0.5);

    constexpr double a1 = 0.5;
    constexpr double a2 = 5.0 /     44.0;
    constexpr double a3 = 1.0 /     66.0;
    constexpr double a4 = 1.0 /    792.0;
    constexpr double a5 = 1.0 /  15840.0;
    constexpr double a6 = 1.0 / 665280.0;

    double x2  = a_x * a_x;
    double x3  = a_x * x2;
    double x4  = x2  * x2;
    double x5  = x2  * x3;
    double x6  = x3  * x3;

    double t1  = a1  * a_x;
    double t2  = a2  * x2;
    double t3  = a3  * x3;
    double t4  = a4  * x4;
    double t5  = a5  * x5;
    double t6  = a6  * x6;

    double se  = 1.0 + t2  + t4 + t6;
    double so  = t1  + t3  + t5;
    double res = (se + so) / (se - so);
    assert(res > 0.0);
    return res;
  }

  //-------------------------------------------------------------------------//
  // ExpPade<long double>: 7th-Order Pade Approximant:                       //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the 7th-order approximant is ~1e-20:
  //
  template<>
  constexpr long double ExpPade<long double>(long double a_x)
  {
    assert(std::fabs(a_x) < 0.5L);

    constexpr long double a1 = 0.5L;
    constexpr long double a2 = 3.0L /       26.0L;
    constexpr long double a3 = 5.0L /      312.0L;
    constexpr long double a4 = 5.0L /     3432.0L;
    constexpr long double a5 = 1.0L /    11440.0L;
    constexpr long double a6 = 1.0L /   308880.0L;
    constexpr long double a7 = 1.0L / 17297280.0L;

    long double x2   = a_x * a_x;
    long double x3   = x2  * a_x;
    long double x4   = x2  * x2;
    long double x5   = x4  * a_x;
    long double x6   = x3  * x3;
    long double x7   = x6  * a_x;

    long double t1   = a1  * a_x;
    long double t2   = a2  * x2;
    long double t3   = a3  * x3;
    long double t4   = a4  * x4;
    long double t5   = a5  * x5;
    long double t6   = a6  * x6;
    long double t7   = a7  * x7;

    long double se   = t6  + t4  + t2 + 1.0L;
    long double so   = t7  + t5  + t3 + t1;
    long double res  = (se + so) / (se - so);
    assert     (res > 0.0L);
    return res;
  }

  //=========================================================================//
  // "LogPade":                                                              //
  //=========================================================================//
  // Pade approximants for log(x), 1/2 <= x < 1, centered at x = 3/4.
  // Altermatively, we could use approximants centered at x=1, which would have
  // "simpler" coeffs, but they would require higher orders for the same preci-
  // sion, so would offer no real advantages, since  in both cases the approxi-
  // mants are "dense":
  //
  template<typename F>
  F LogPade(F a_x);

  //-------------------------------------------------------------------------//
  // LogPade<float>: 4th-Order Pade Approximant:                             //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-23) =~ 1.2e-7,
  // the abs precision of the 4th-order approximant is ~6e-9 (yet, the 3rd-order
  // one would be insufficient):
  //
  template<>
  constexpr float LogPade<float>(float a_x)
  {
    assert(0.5f <= a_x && a_x < 1.0f);

    constexpr float  b4 =  1536.0f;
    constexpr float  b3 = 18432.0f;
    constexpr float  b2 = 31104.0f;
    constexpr float  b1 = 10368.0f;
    constexpr float  b0 =   486.0f;

    constexpr double A  = -0.28768207245178093;           // log(3/4)
    constexpr float  a4 = float(double(b4) * A +  6400.0);
    constexpr float  a3 = float(double(b3) * A + 30720.0);
    constexpr float  a2 = float(double(b2) * A);
    constexpr float  a1 = float(double(b1) * A - 17280.0);
    constexpr float  a0 = float(double(b0) * A -  2025.0);

    float res =
      ((((a4 * a_x + a3) * a_x + a2) * a_x + a1) * a_x + a0) /
      ((((b4 * a_x + b3) * a_x + b2) * a_x + b1) * a_x + b0);
    assert(res < 0.0f);
    return res;
  }

  //-------------------------------------------------------------------------//
  // LogPade<double>: 8th-Order Pade Approximant:                            //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-52) =~ 2.2e-16,
  // the abs precision of the 8th-order approximation is ~7e-17:
  //
  template<>
  constexpr double LogPade<double>(double a_x)
  {
    assert(0.5 <= a_x && a_x < 1.0);

    constexpr double b8 =     9175040.0;
    constexpr double b7 =   440401920.0;
    constexpr double b6 =  4046192640.0;
    constexpr double b5 = 12138577920.0;
    constexpr double b4 = 14224896000.0;
    constexpr double b3 =  6827950080.0;
    constexpr double b2 =  1280240640.0;
    constexpr double b1 =    78382080.0;
    constexpr double b0 =      918540.0;

    constexpr double A  =  -0.28768207245178093;          // log(3/4)
    constexpr double a8 =  b8 * A +    49872896.0;
    constexpr double a7 =  b7 * A +  1402994688.0;
    constexpr double a6 =  b6 * A +  7687766016.0;
    constexpr double a5 =  b5 * A + 10924720128.0;
    constexpr double a4 =  b4 * A;
    constexpr double a3 =  b3 * A -  6145155072.0;
    constexpr double a2 =  b2 * A -  2432457216.0;
    constexpr double a1 =  b1 * A -   249702912.0;
    constexpr double a0 =  b0 * A -     4992921.0;

    double res =
      ((((((((a8 * a_x + a7) * a_x + a6) * a_x + a5) * a_x + a4) * a_x + a3)
                 * a_x + a2) * a_x + a1) * a_x + a0) /
      ((((((((b8 * a_x + b7) * a_x + b6) * a_x + b5) * a_x + b4) * a_x + b3)
                 * a_x + b2) * a_x + b1) * a_x + b0);
    assert(res < 0.0L);
    return res;
  }

  //-------------------------------------------------------------------------//
  // LogPade<long double>: 10th-Order Pade Approximant:                      //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the 10th-order approximant is ~8e-21:
  //
  template<>
  constexpr long double LogPade<long double>(long double a_x)
  {
    assert(0.5L  <= a_x && a_x < 1.0L);

    constexpr long double b10 =     1321205760.0L;
    constexpr long double b9  =    99090432000.0L;
    constexpr long double b8  =  1504935936000.0L;
    constexpr long double b7  =  8026324992000.0L;
    constexpr long double b6  = 18435465216000.0L;
    constexpr long double b5  = 19910302433280.0L;
    constexpr long double b4  = 10369949184000.0L;
    constexpr long double b3  =  2539579392000.0L;
    constexpr long double b2  =   267846264000.0L;
    constexpr long double b1  =     9920232000.0L;
    constexpr long double b0  =       74401740.0L;

    constexpr long double A   = -0.287682072451780927439L;   // log(3/4)
    constexpr long double a10 = b10 * A +     7739539456.0L;
    constexpr long double a9  = b9  * A +   362466508800.0L;
    constexpr long double a8  = b8  * A +  3665593958400.0L;
    constexpr long double a7  = b7  * A + 12192369868800.0L;
    constexpr long double a6  = b6  * A + 13519341158400.0L;
    constexpr long double a5  = b5  * A;
    constexpr long double a4  = b4  * A -  7604629401600.0L;
    constexpr long double a3  = b3  * A -  3857742028800.0L;
    constexpr long double a2  = b2  * A -   652396971600.0L;
    constexpr long double a1  = b1  * A -    36287578800.0L;
    constexpr long double a0  = b0  * A -      435840669.0L;

    long double res =
      ((((((((((a10 * a_x + a9) * a_x + a8) * a_x + a7) * a_x + a6) * a_x + a5)
                    * a_x + a4) * a_x + a3) * a_x + a2) * a_x + a1) * a_x + a0)
      /
      ((((((((((b10 * a_x + b9) * a_x + b8) * a_x + b7) * a_x + b6) * a_x + b5)
                    * a_x + b4) * a_x + b3) * a_x + b2) * a_x + b1) * a_x + b0);
    assert(res < 0.0L);
    return res;
  }

  //=========================================================================//
  // "CosPade":                                                              //
  //=========================================================================//
  // Assumes 0 <= x <= Pi/4. We have a choice of 2 kinds of approximants here:
  // (a) those centered at x=0,   containing polynomials of even degrees only;
  // (b) those centered at x=Pi/8 with "dense" polynomials  and  quite complex
  //     coeffs, but of slighly lower degree;
  // from the computational efficiency point of view, method (a) is better!
  //
  template<typename F>
  F CosPade(F a_x);

  //-------------------------------------------------------------------------//
  // CosPade<float>: 4th-Order Pade approximant:                             //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-23) =~ 1.2e-7,
  // the abs precision of the 4th-order approximant is ~4e-8:
  //
  template<>
  constexpr float CosPade<float>(float a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<float> + 10.0f * Eps<float>);

    constexpr float a4 = float( 313.0 / 15120.0);
    constexpr float a2 = float(-115.0 /   252.0);
    constexpr float b4 = float(  13.0 / 15120.0);
    constexpr float b2 = float(  11.0 /   252.0);
    float           x2 = a_x * a_x;
    return
      ((a4 * x2 + a2) * x2 + 1.0f) /
      ((b4 * x2 + b2) * x2 + 1.0f);
  }

  //-------------------------------------------------------------------------//
  // CosPade<double>: 8th-Order Pade approximant:                            //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-52) =~ 2.2e-16,
  // the abs precision of the 8th-order approximation is ~2e-18 (however, the
  // 6th-order one would be insufficient):
  //
  template<>
  constexpr double CosPade(double a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<double> + 100.0 * Eps<double>);

    constexpr double a8 =   80737373.0 / 23594700729600.0;
    constexpr double a6 = -  7696415.0 /    13108167072.0;
    constexpr double a4 =    4375409.0 /      141863280.0;
    constexpr double a2 = -   260735.0 /         545628.0;

    constexpr double b8 =      11321.0 /  1814976979200.0;
    constexpr double b6 =     109247.0 /    65540835360.0;
    constexpr double b4 =      34709.0 /      141863280.0;
    constexpr double b2 =      12079.0 /         545628.0;
    double           x2 = a_x * a_x;
    return
      ((((a8 * x2 + a6) * x2 + a4) * x2 + a2) * x2 + 1.0) /
      ((((b8 * x2 + b6) * x2 + b4) * x2 + b2) * x2 + 1.0);
  }

  //-------------------------------------------------------------------------//
  // CosPade<long double>: (10,8)th-Order Pade approximant:                  //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the (10,8)th-order approximant is 1.5e-21:
  //
  template<>
  constexpr long double CosPade(long double a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<long double> + 100.0L * Eps<long double>);

    constexpr long double a10 = -  19934927.0L / 865535480409600.0L;
    constexpr long double a8  =   106649201.0L /  15274155536640.0L;
    constexpr long double a6  = - 250412863.0L /    327854348640.0L;
    constexpr long double a4  =     1341515.0L /        40031056.0L;
    constexpr long double a2  = -   7257019.0L /        15011646.0L;

    constexpr long double b8  =      391693.0L / 259660644122880.0L;
    constexpr long double b6  =      101797.0L /    163927174320.0L;
    constexpr long double b4  =        5293.0L /        40031056.0L;
    constexpr long double b2  =      124402.0L /         7505823.0L;
    long double           x2  = a_x * a_x;
    return
      (((((a10 * x2 + a8) * x2 + a6) * x2 + a4) * x2 + a2) * x2 + 1.0L) /
       (((( b8 * x2 + b6) * x2 + b4) * x2 + b2) * x2 + 1.0L);
  }

  //=========================================================================//
  // "SinPade":                                                              //
  //=========================================================================//
  // The arg is assumed to be in [0 .. Pi/4], centered at x=0 (same considera-
  // tions apply as for "CosPade"):
  //
  template<typename F>
  F SinPade(F a_x);

  //-------------------------------------------------------------------------//
  // SinPade<float>: (5,4)th-Order Pade approximant:                         //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-23) =~ 1.2e-7,
  // the abs precision of the (5,4)th-order approximant is ~1.6e-9:
  //
  template<>
  constexpr float SinPade<float>(float a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<float> + float(10.0 * Eps<float>));
    constexpr float a4 = float( 551.0 / 166320.0);
    constexpr float a2 = float(- 53.0 /    396.0);
    constexpr float b4 = float(   5.0 /  11088.0);
    constexpr float b2 = float(  13.0 /    396.0);
    float           x2 = a_x * a_x;
    return
      a_x * ((a4 * x2 + a2) * x2 + 1.0f) /
            ((b4 * x2 + b2) * x2 + 1.0f);
  }

  //-------------------------------------------------------------------------//
  // SinPade<double>: (7,8)th-Order Pade approximant:                        //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-52) =~ 2.2e-16,
  // the abs precision of the (7,8)th-order approximation is ~7e-17:
  //
  template<>
  constexpr double SinPade<double>(double a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<double> + 100.0 * Eps<double>);

    constexpr double a6 = - 62077121.0 /   1727021696400.0;
    constexpr double a4 =    9713777.0 /      2242885320.0;
    constexpr double a2 = -  2020961.0 /        14377470.0;

    constexpr double b8 =    1768969.0 / 124345562140800.0;
    constexpr double b6 =      36317.0 /     12335869260.0;
    constexpr double b4 =      26015.0 /        74762844.0;
    constexpr double b2 =     187642.0 /         7188735.0;
    double           x2 = a_x * a_x;
    return
      a_x *  (((a6 * x2 + a4) * x2 + a2) * x2 + 1.0)          /
            ((((b8 * x2 + b6) * x2 + b4) * x2 + b2) * x2 + 1.0);
  }

  //-------------------------------------------------------------------------//
  // SinPade<long double>: (10,8)th-Order Pade approximant:                  //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the (10,8)th-order approximant is about the same:
  //
  template<>
  constexpr long double SinPade<long double>(long double a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<long double> + 100.0L * Eps<long double>);

    constexpr long double a8 =  4585922449.0L / 15605159573203200.0L;
    constexpr long double a6 = - 269197963.0L /     3940696861920.0L;
    constexpr long double a4 =    38518909.0L /        7217393520.0L;
    constexpr long double a2 = -  53272705.0L /         360869676.0L;

    constexpr long double b8 =     1029037.0L /   346781323848960.0L;
    constexpr long double b6 =      560401.0L /      562956694560.0L;
    constexpr long double b4 =     1281433.0L /        7217393520.0L;
    constexpr long double b2 =     2290747.0L /         120289892.0L;
    long double           x2 = a_x * a_x;
    return
      a_x * ((((a8 * x2 + a6) * x2 + a4) * x2 + a2) * x2 + 1.0L) /
            ((((b8 * x2 + b6) * x2 + b4) * x2 + b2) * x2 + 1.0L);
  }

  //=========================================================================//
  // "Exp", "Log", "Pow", "Cos", "Sin" for arbitrary real args:              //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "Exp" for an arbitrary real arg:                                        //
  //-------------------------------------------------------------------------//
  template<bool ForceApprox=true, typename F>
  constexpr F Exp(F a_x)
  {
    // Complex "Exp" has a separate specific implementation:
    static_assert(!IsComplex<F>);

    // Special Cases:
    if (std::isnan(a_x))
      return NaN<F>;
    if (std::isinf(a_x))
      return (a_x > 0) ? Inf<F> : F(0);

    // Generic Case:
#   ifndef __clang__
    if constexpr(ForceApprox)
    {
#   endif
      // First, change the base from E to 2, to make it easier to represent the
      // result in the IEEE 754 format: exp(x) = 2^y, where y = x * Log2E:
      F y = a_x * Log2E<F>;

      // Get the integral and fractional part of "y":   y = intgY + fracY:
      F intgY = NaN<F>;
      F fracY = std::modf (y, &intgY);
      assert(std::isfinite(intgY) && std::fabs(fracY) < F(1.0));

      // Check if "intgY" is too large in absolute value:
      if (intgY > F(INT_MAX))
        return Inf<F>;
      if (intgY < F(INT_MIN))
        return F(0);

      // Generic Case:
      int n = int(intgY);

      // Make "fracY" within [-0.5 .. 0.5] for better convergence:
      if (fracY < F(-0.5))
      {
        fracY += F(1.0);
        --n;
      }
      else
      if (fracY > F(0.5))
      {
        fracY -= F(1.0);
        ++n;
      }
      assert(std::fabs(fracY) <= F(0.5));

      // For the fractional part, get back to Base-E; it can only become smaller
      // in abs value as a result:
      F f = fracY * Ln2<F>;
      assert(std::fabs(f) <  F(0.5));

      // Use the Pade approximant for exp(f) around f=0:
      F res = ExpPade<F>(f);

      // Multiply the result by 2^n (can still get 0 or +oo at this moment):
      return std::ldexp(res, n);

#   ifndef __clang__
    }
    else
      // GCC, and NOT using the Pade approximant:
      return std::exp(a_x);
#   endif
  }

  //-------------------------------------------------------------------------//
  // "Log" (Natural Logarithm) for an arbitrary real arg:                    //
  //-------------------------------------------------------------------------//
  template<bool ForceApprox=true, typename F>
  constexpr F Log(F a_x)
  {
    // Complex "Log" has a separate specific implementation:
    static_assert(!IsComplex<F>);

    // Special Cases:
    if (std::isnan(a_x) || a_x < 0)
      return NaN<F>;

    if (a_x == 0)
      return -Inf<F>;

    if (std::isinf(a_x))
    {
      assert(a_x > 0);
      return  Inf<F>;
    }

    // Generic Case:
#   ifndef __clang__
    if constexpr(ForceApprox)
    {
#   endif
      // Get the Base-2 exponent and the normalised fractional part of "a_x":
      assert(std::isfinite(a_x) && a_x > 0);
      int e2X   = INT_MIN;
      F   fracX = std::frexp(a_x,  &e2X);
      assert(0.5 <= fracX && fracX < F(1) && e2X != INT_MIN);

      // Then log(x) = log(2^e2X  * fracX) = e2X * log(2) + log(fracX),
      // so expand     log(fracX) in the vicinity of fracX = 3/4:
      F logFX = LogPade<F>(fracX);
      assert(logFX < 0);

      // From "e2X" (which is actually the main part of the Base-2 log of
      // "a_x"), get back to the natural log:
      return F(e2X) * Ln2<F> + logFX;

#   ifndef __clang__
    }
    else
      // GCC, and NOT using the Pade approximant:
      return std::log(a_x);
#   endif
  }

  //-------------------------------------------------------------------------//
  // Generic "Pow" for an arbitrary real arg:                                //
  //-------------------------------------------------------------------------//
  template<bool ForceApprox=true, typename F>
  constexpr F Pow(F a_x, F a_y)
  {
    // Complex "Pow" is not implemented yet. Also, for this generic implementa-
    // tion,   "a_x" is required to be positive:
    static_assert(!IsComplex<F>);
    assert(a_x > 0);
    return Exp<F>(a_y * Log<F>(a_x));
  }

  //-------------------------------------------------------------------------//
  // "Cos" for an arbitrary real arg:                                        //
  //-------------------------------------------------------------------------//
  template<bool ForceApprox=true, typename F>
  constexpr F Cos(F a_x)
  {
    // Complex "Cos" has a separate specific implementation:
    static_assert(!IsComplex<F>);

    if (!std::isfinite(a_x))
      return NaN<F>;

    // Use the Pade approximant if the compiler is CLang   (because in that
    // case, "std::cos" is not "constexpr"), OR if explictly asked to do so.
    // But the Pade approximant is not applicable to complex args:
#   ifndef __clang__
    if constexpr(ForceApprox)
    {
#   endif
      // First, normalise the arg to the interval [0; +00):
      // NB: "fabs" and "fmod" are "constexpr" functions in C++ >= 23:
      a_x = std::fabs(a_x);

      // Then normalise it to the interval [0..2*Pi):
      a_x = std::fmod(a_x, TwoPi<F>);

      // Then to the interval [0..Pi]:
      bool chSgn = (a_x > Pi<F>);
      if  (chSgn)
        a_x -= Pi<F>;

      // Then to the interval [0..Pi/2]:
      if (a_x > Pi_2<F>)
      {
        a_x   = Pi<F> - a_x;
        chSgn = !chSgn;
      }
      // Finally, to the interval [0..Pi/4], and compute the function:
      F res = (a_x <= Pi_4<F>) ? CosPade<F>(a_x) : SinPade<F>(Pi_2<F> - a_x);

      // Don't forget the sign:
      if (chSgn)
        res = - res;
      assert(std::fabs(res) < F(1.0) + Eps<F>);
      return res;
#   ifndef __clang__
    }
    else
      // GCC, and NOT using the Pade approximant:
      return std::cos(a_x);
#   endif
  }

  //-------------------------------------------------------------------------//
  // "Sin" for an arbitrary real arg:                                        //
  //-------------------------------------------------------------------------//
  template<bool ForceApprox=true, typename F>
  constexpr F Sin(F a_x)
  {
    // Complex "Sin" has a separate specific implementation:
    static_assert(!IsComplex<F>);

    if (!std::isfinite(a_x))
      return NaN<F>;

    // Use the Pade approximant if the compiler is CLang   (because in that
    // case, "std::sin" is not "constexpr"), OR if explictly asked to do so.
    // But the Pade approximant is not applicable to complex args:
#   ifndef __clang__
    if constexpr(ForceApprox)
    {
#   endif
      // First, normalise the arg to the interval [0; +00):
      // NB: "fabs" and "fmod" are "constexpr" functions in C++ >= 23:
      bool  chSgn =  (a_x < 0);
      a_x = std::fabs(a_x);

      // Then normalise it to the interval [0..2*Pi):
      a_x = std::fmod(a_x, TwoPi<F>);

      // Then to the interval [0..Pi]:
      if (a_x > Pi<F>)
      {
        a_x  -= Pi<F>;
        chSgn = !chSgn;
      }
      // Then to the interval [0..Pi/2] (same sign):
      if (a_x > Pi_2<F>)
        a_x = Pi<F> - a_x;

      // Finally, to the interval [0..Pi/4], and compute the function:
      F res = (a_x <= Pi_4<F>) ? SinPade<F>(a_x) : CosPade<F>(Pi_2<F> - a_x);

      // Don't forget the sign:
      if (chSgn)
        res = - res;
      assert(std::fabs(res) < F(1.0) + Eps<F>);
      return res;
#   ifndef __clang__
    }
    else
      // GCC, and NOT using the Pade approximant:
      return std::sin(a_x);
#   endif
  }

  //=========================================================================//
  // Square and Cubic Roots:                                                 //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "SqRt":                                                                 //
  //-------------------------------------------------------------------------//
  // FIXME: This is a placeholder: this implementation is not "constexpr" in
  // CLang:
  template<typename F>
  constexpr F SqRt (F a_x)
  {
    if (!(a_x >= 0))
      return NaN<F>;
    return std::sqrt (a_x);
  }

  //-------------------------------------------------------------------------//
  // "CbRt":                                                                 //
  //-------------------------------------------------------------------------//
  // FIXME: This is a placeholder as well:
  //
  template<typename F>
  constexpr F CbRt (F a_x)
  {
    return std::cbrt (a_x);
  }

  //=========================================================================//
  // Complex Functions:                                                      //
  //=========================================================================//
  // The templates below specifically match "std::complex" args.
  // The args are passed by copy -- for relatively small "std::complex" objs,
  // this is most efficient:
  //
  //-------------------------------------------------------------------------//
  // Complex "Exp":                                                          //
  //-------------------------------------------------------------------------//
  template<typename T>
  constexpr std::complex<T> Exp(std::complex<T> a_z)
  {
    // exp(x+I*y) = exp(x)*(cos(y)+I*sin(y)):
    T x     = std::real(a_z);
    T y     = std::imag(a_z);
    T r     = Exp(x);
    return std::complex<T>(r * Cos(y), r * Sin(y));
  }

  //-------------------------------------------------------------------------//
  // Complex "Cos":                                                          //
  //-------------------------------------------------------------------------//
  template<typename T>
  constexpr std::complex<T> Cos(std::complex<T> a_z)
  {
    // cos(x+I*y) = cos(x)*ch(y) - I*sin(x)*sh(y):
    T x     =   std::real(a_z);
    T y     =   std::imag(a_z);
    T u     =   Exp(y);
    T u1    =   T(1)     / u;
    T chY   =   (u + u1) / T(2);
    T shY   =   (u - u1) / T(2);
    T reCos =   Cos(x) * chY;
    T imCos = - Sin(x) * shY;
    return std::complex<T>(reCos, imCos);
  }

  //-------------------------------------------------------------------------//
  // Complex "Sin":                                                          //
  //-------------------------------------------------------------------------//
  template<typename T>
  constexpr std::complex<T> Sin(std::complex<T> a_z)
  {
    // sin(x+I*y) = sin(x)*ch(y) + I*cos(x)*sh(y):
    T x     = std::real(a_z);
    T y     = std::imag(a_z);
    T u     = Exp(y);
    T u1    = T(1)     / u;
    T chY   = (u + u1) / T(2);
    T shY   = (u - u1) / T(2);
    T reSin = Sin(x) * chY;
    T imSin = Cos(x) * shY;
    return std::complex<T>(reSin, imSin);
  }

  //-------------------------------------------------------------------------//
  // Complex "Cos" and "Sin" together (for efficiency):                      //
  //-------------------------------------------------------------------------//
  // Returns (cos(z), sin(z)):
  //
  template<typename T>
  constexpr std::pair<std::complex<T>, std::complex<T>>
    CosSin (std::complex<T> a_z)
  {
    T x     =   std::real(a_z);
    T y     =   std::imag(a_z);
    T u     =   Exp(y);
    T u1    =   T(1)     / u;
    T chY   =   (u + u1) / T(2);
    T shY   =   (u - u1) / T(2);
    T cosX  =   Cos(x);
    T sinX  =   Sin(x);
    T reCos =   cosX * chY;
    T imCos = - sinX * shY;
    T reSin =   sinX * chY;
    T imSin =   cosX * shY;
    std::complex<T>  cosZ(reCos, imCos);
    std::complex<T>  sinZ(reSin, imSin);
    return std::make_pair(cosZ,  sinZ) ;
  }

  //-------------------------------------------------------------------------//
  // Complex "Pow":                                                          //
  //-------------------------------------------------------------------------//
  // FIXME: This implementation is NOT actually "constexpr" in CLang:
  //
  template<typename T>
  constexpr  std::complex<T> Pow(std::complex<T> a_z, T a_p)
    { return std::pow(a_z, a_p); }

  //-------------------------------------------------------------------------//
  // Complex "SqRt":                                                         //
  //-------------------------------------------------------------------------//
  // FIXME: Also a placeholder as yet:
  //
  template<typename T>
  constexpr  std::complex<T> SqRt(std::complex<T> a_z)
    { return std::sqrt(a_z); }

  //-------------------------------------------------------------------------//
  // Complex "CbRt":                                                         //
  //-------------------------------------------------------------------------//
  // FIXME: Also as placeholder as yet:
  //
  template<typename T>
  constexpr  std::complex<T> CbRt(std::complex<T> a_z)
    { return std::pow(a_z, T(1)/T(3)); }
}
