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
#include <type_traits>
#include <iostream>

namespace DimTypes::Bits::CEMaths
{
  //=========================================================================//
  // Elementary Functions, Utils:                                            //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "IsComplex":                                                            //
  //-------------------------------------------------------------------------//
  template<typename T> inline constexpr bool IsComplex = false;
  template<typename T> inline constexpr bool IsComplex<std::complex<T>> = true;

  //-------------------------------------------------------------------------//
  // "Abs":                                                                  //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F Abs(F a_x)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    return std::abs(a_x);
  }

  //-------------------------------------------------------------------------//
  // "Floor":                                                                //
  //-------------------------------------------------------------------------//
  // XXX: Not really "constexpr" in CLang <= 19, and no emulation is currently
  // provided (TODO: not very difficult to implement):
  //
  template<typename F>
  constexpr F Floor(F a_x)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    return std::floor(a_x);
  }

  //-------------------------------------------------------------------------//
  // "Ceil":                                                                 //
  //-------------------------------------------------------------------------//
  // XXX: Not really "constexpr" in CLang <= 19, and no emulation is currently
  // provided (TODO: not very difficult to implement):
  //
  template<typename F>
  constexpr F Ceil(F a_x)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    return std::ceil(a_x);
  }

  //-------------------------------------------------------------------------//
  // "Round":                                                                //
  //-------------------------------------------------------------------------//
  // XXX: Not really "constexpr" in CLang <= 19, and no emulation is currently
  // provided (TODO: not very difficult to implement):
  //
  template<typename F>
  constexpr F Round(F a_x)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    return std::round(a_x);
  }

  //-------------------------------------------------------------------------//
  // "FMod":                                                                 //
  //-------------------------------------------------------------------------//
  // Again, in CLang <= 19, "std::fmod" is not "constexpr", so our own implem-
  // entation is needed:
  //
  template<typename F>
  constexpr F FMod (F a_x, F a_y)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    // For the avoidance of doubt, we don't support negative args:
    assert(a_x >= 0.0 && a_y > 0.0);
    F q  = a_x / a_y;
    assert(q >= 0.0);
    // Convert "q" to "unsigned long" -- it is always rounded down:
    using  ULong = unsigned long;
    ULong  u = ULong(q);
    F      r = std::max(a_x - F(u) * a_y, F(0.0));
    return r;
# else
    return std::fmod(a_x, a_y);
# endif
  }

  // NB: "std::isinf", "std::isnan", "std::isfinite" are ALWAYS "constexpr" in
  // compilers supporting C++ >= 23 (even in CLang <= 19 which is not entirely
  // standard-complient), so there is no need to wrap them into our own impls...

  //-------------------------------------------------------------------------//
  // "Sqr", "Cube":                                                          //
  //-------------------------------------------------------------------------//
  template<typename F> constexpr F Sqr (F a_x) { return a_x * a_x;       }
  template<typename F> constexpr F Cube(F a_x) { return a_x * a_x * a_x; }

  //=========================================================================//
  // Consts:                                                                 //
  //=========================================================================//
  // Infinity, NaN:
  template<typename F> inline constexpr F Inf   = F(INFINITY);     // +oo
  template<typename F> inline constexpr F NaN   = F(NAN);          // NaN

  // Eps:
  template<typename F> inline constexpr F Eps;
  template<> inline constexpr float       Eps<float>       = FLT_EPSILON;
  template<> inline constexpr double      Eps<double>      = DBL_EPSILON;
  template<> inline constexpr long double Eps<long double> = LDBL_EPSILON;

  // Pi:
  template<typename F> inline constexpr F Pi      = F(M_PIl);        // Pi
  template<typename F> inline constexpr F TwoPi   = F(2.0L * M_PIl); // 2*Pi
  template<typename F> inline constexpr F Pi_2    = F(M_PI_2l);      // Pi/2
  template<typename F> inline constexpr F Pi_4    = F(M_PI_4l);      // Pi/4

  // "Technical" Consts:
  template<typename F> inline constexpr F SqRt2   = F(M_SQRT2l);     // SqRt(2)
  template<typename F> inline constexpr F SqRt1_2 = F(M_SQRT1_2l);   // .../2
  template<typename F> inline constexpr F Ln2     = F(M_LN2l);       // ln(2)
  template<typename F> inline constexpr F Log2E   = F(M_LOG2El);     // log2(e)

  template<typename F> inline constexpr F SqRt3;                     // SqRt(3)
  template<typename F> inline constexpr F CbRt2;                     // CbRt(2)
  template<typename F> inline constexpr F CbRt4;                     // CbRt(4)
  template<typename F> inline constexpr F CbRt48;                    // CbRt(48)
  template<typename F> inline constexpr F Ln3_4;                     // Ln(3/4),

  // XXX: Suppress bogus warnings generated by CLang <= 19 on the following
  // constant specialisations:
  template<> inline constexpr long double SqRt3<long double>  =
    1.73205080756887729353L;
  template<> inline constexpr double      SqRt3<double>       =
    double(SqRt3<long  double>);
  template<> inline constexpr float       SqRt3<float>        =
    float (SqRt3<long  double>);

  template<> inline constexpr long double CbRt2<long double>  =
    1.25992104989487316477L;
  template<> inline constexpr double      CbRt2<double>       =
    double(CbRt2<long double>);
  template<> inline constexpr float       CbRt2<float>        =
    float (CbRt2<long double>);

  template<> inline constexpr long double CbRt4<long double>  =
    1.58740105196819947475L;
  template<> inline constexpr double      CbRt4<double>       =
    double(CbRt4<long double>);
  template<> inline constexpr float       CbRt4<float>        =
    float (CbRt4<long double>);

  template<> inline constexpr long double CbRt48<long double> =
    3.63424118566427931778L;
  template<> inline constexpr double      CbRt48<double>      =
    double(CbRt48<long double>);
  // CbRt48<float> is proibably not required

  template<> inline constexpr long double Ln3_4<long double>  =
    -0.287682072451780927439L;
  template<> inline constexpr double      Ln3_4<double>       =
    double(Ln3_4<long  double>);
  // Ln3_4<float> is proibably not required

  //-------------------------------------------------------------------------//
  // "ApproxEqual":                                                          //
  //-------------------------------------------------------------------------//
  template<typename F>
  inline constexpr  F DefaultTol = Eps<F> * F(100.0);

  template<typename F>
  constexpr bool ApproxEqual(F a_x, F a_y, F a_tol = DefaultTol<F>)
  {
    assert(a_tol >= F(0.0));
    F err =
      Abs
      (
        Abs(a_y) < F(1.0)
        ? a_x - a_y
        : a_x / a_y - F(1.0)
      );
    return err < a_tol;
  }

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
  constexpr F ExpPade(F a_x);

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
    assert(Abs(a_x) < 0.5f);

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
    assert(Abs(a_x) < 0.5);

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
    assert(Abs(a_x) < 0.5L);

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
  constexpr F LogPade(F a_x);

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

    constexpr float  a4 = float(double(b4) * Ln3_4<double> +  6400.0);
    constexpr float  a3 = float(double(b3) * Ln3_4<double> + 30720.0);
    constexpr float  a2 = float(double(b2) * Ln3_4<double>);
    constexpr float  a1 = float(double(b1) * Ln3_4<double> - 17280.0);
    constexpr float  a0 = float(double(b0) * Ln3_4<double> -  2025.0);

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

    constexpr double a8 =  b8 * Ln3_4<double> +    49872896.0;
    constexpr double a7 =  b7 * Ln3_4<double> +  1402994688.0;
    constexpr double a6 =  b6 * Ln3_4<double> +  7687766016.0;
    constexpr double a5 =  b5 * Ln3_4<double> + 10924720128.0;
    constexpr double a4 =  b4 * Ln3_4<double>;
    constexpr double a3 =  b3 * Ln3_4<double> -  6145155072.0;
    constexpr double a2 =  b2 * Ln3_4<double> -  2432457216.0;
    constexpr double a1 =  b1 * Ln3_4<double> -   249702912.0;
    constexpr double a0 =  b0 * Ln3_4<double> -     4992921.0;

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

    constexpr long double a10 = b10 * Ln3_4<long double> +     7739539456.0L;
    constexpr long double a9  = b9  * Ln3_4<long double> +   362466508800.0L;
    constexpr long double a8  = b8  * Ln3_4<long double> +  3665593958400.0L;
    constexpr long double a7  = b7  * Ln3_4<long double> + 12192369868800.0L;
    constexpr long double a6  = b6  * Ln3_4<long double> + 13519341158400.0L;
    constexpr long double a5  = b5  * Ln3_4<long double>;
    constexpr long double a4  = b4  * Ln3_4<long double> -  7604629401600.0L;
    constexpr long double a3  = b3  * Ln3_4<long double> -  3857742028800.0L;
    constexpr long double a2  = b2  * Ln3_4<long double> -   652396971600.0L;
    constexpr long double a1  = b1  * Ln3_4<long double> -    36287578800.0L;
    constexpr long double a0  = b0  * Ln3_4<long double> -      435840669.0L;

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
  constexpr F CosPade(F a_x);

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
  constexpr F SinPade(F a_x);

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
  // SinPade<long double>: (9,8)th-Order Pade approximant:                   //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the (9,8)th-order approximant is about the same:
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
  // "ATanPade":                                                             //
  //=========================================================================//
  // Pade approximants for atan(x), |x| <= 1, centered at x = 0;
  // Compared to other elementary functions, atan(x) requires relatively
  // lengthy Pade approximants to achieve the required accuracy:
  // XXX: A rational Pade approximation may not be the most efficient choice in
  // this case; Chebyshev approximation constructed using the Remez algorithm
  // might be better. However, it may not preserbe the exact value and derivat-
  // ives at x=0. Anoher possibiliy is to use a muli-point atan approximation:
  //
  template<typename F>
  constexpr F ATanPade(F a_x);

  //-------------------------------------------------------------------------//
  // ATanPade<float>: (11,8)th-Order Pade Approximant:                       //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-23) =~ 1.2e-7,
  // the abs precision of the (11,8)th-order approximant is ~3.7e-8:
  //
  template<>
  constexpr float ATanPade<float>(float a_x)
  {
    assert(-1.0f <= a_x && a_x <= 1.0f);

    constexpr float a11 = float(- 16384.0 / 72747675.0);
    constexpr float  a9 = float(  16384.0 /  1322685.0);
    constexpr float  a7 = float(   3159.0 /    11305.0);
    constexpr float  a5 = float(   6139.0 /     4845.0);
    constexpr float  a3 = float(    113.0 /       57.0);

    constexpr float  b8 = float(    231.0 /     4199.0);
    constexpr float  b6 = float(    924.0 /     1615.0);
    constexpr float  b4 = float(    594.0 /      323.0);
    constexpr float  b2 = float(     44.0 /       19.0);

    float            x2 = a_x * a_x;
    return
      (((((a11 * x2 + a9) * x2 + a7) * x2 + a5) * x2 + a3) * x2 + 1.0f) * a_x /
       ((((b8  * x2 + b6) * x2 + b4) * x2 + b2) * x2 + 1.0f);
  }

  //-------------------------------------------------------------------------//
  // ATanPade<double>: (21,20)th-Order Pade Approximant:                     //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-52) =~ 2.2e-16,
  // the abs precision of the (21,20)th-order approximant is ~1.1e-16:
  //
  template<>
  constexpr double ATanPade<double>(double a_x)
  {
    assert(-1.0 <= a_x && a_x <= 1.0);

    constexpr double a21 = 68719476736.0 / 65261681526586545.0;
    constexpr double a19 =   562144147.0 /     2456679146493.0;
    constexpr double a17 =   350944637.0 /       43099634149.0;
    constexpr double a15 =   326695412.0 /        2925314535.0;
    constexpr double a13 =      140612.0 /            181753.0;
    constexpr double a11 =      153386.0 /             50061.0;
    constexpr double  a9 =      100454.0 /             13653.0;
    constexpr double  a7 =     1504228.0 /            138047.0;
    constexpr double  a5 =       25908.0 /              2665.0;
    constexpr double  a3 =         589.0 /               123.0;

    constexpr double b20 =        2261.0 /         156835045.0;
    constexpr double b18 =        4522.0 /           4091349.0;
    constexpr double b16 =       33915.0 /           1363783.0;
    constexpr double b14 =      348840.0 /           1363783.0;
    constexpr double b12 =       67830.0 /             47027.0;
    constexpr double b10 =       81396.0 /             16687.0;
    constexpr double  b8 =      203490.0 /             19721.0;
    constexpr double  b6 =      271320.0 /             19721.0;
    constexpr double  b4 =        5985.0 /               533.0;
    constexpr double  b2 =         210.0 /                41.0;

    double            x2 = a_x * a_x;
    return
      ((((((((((a21 * x2 + a19) * x2 + a17) * x2 + a15) * x2 + a13) * x2 + a11)
                    * x2 +  a9) * x2 +  a7) * x2 +  a5) * x2 +  a3) * x2 + 1.0)
      * a_x /
      ((((((((((b20 * x2 + b18) * x2 + b16) * x2 + b14) * x2 + b12) * x2 + b10)
                    * x2 +  b8) * x2 +  b6) * x2 +  b4) * x2 +  b2) * x2 + 1.0);
  }

  //-------------------------------------------------------------------------//
  // ATanPade<long double>: (27,24)th-Order Pade Approximant:                //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the (27,24)th-order approximant is ~1.8e-20:
  //
  template<>
  constexpr long double ATanPade<long double>(long double a_x)
  {
    assert(-1.0L <= a_x && a_x <= 1.0L);

    constexpr long double a27 = -17592186044416.0L / 13619948432012945752275.0L;
    constexpr long double a25 =  17592186044416.0L /    38803271886076768525.0L;
    constexpr long double a23 =    567661032929.0L /        7645964903660447.0L;
    constexpr long double a21 =     86486630155.0L /          32170960323957.0L;
    constexpr long double a19 =    456530944307.0L /          10723653441319.0L;
    constexpr long double a17 =     61212140501.0L /            166000827265.0L;
    constexpr long double a15 =      1543545098.0L /               791736855.0L;
    constexpr long double a13 =     27146901662.0L /              4064249189.0L;
    constexpr long double a11 =       116611930.0L /                 7625233.0L;
    constexpr long double  a9 =         8318998.0L /                  352359.0L;
    constexpr long double  a7 =         4750009.0L /                  195755.0L;
    constexpr long double  a5 =           66263.0L /                    4165.0L;
    constexpr long double  a3 =             307.0L /                      51.0L;

    constexpr long double b24 =           98325.0L /             19293438101.0L;
    constexpr long double b22 =          235980.0L /               665290969.0L;
    constexpr long double b20 =         1297890.0L /               150226993.0L;
    constexpr long double b18 =        15863100.0L /               150226993.0L;
    constexpr long double b16 =        16223625.0L /                21460999.0L;
    constexpr long double b14 =        25957800.0L /                 7540351.0L;
    constexpr long double b12 =          865260.0L /                   82861.0L;
    constexpr long double b10 =         5191560.0L /                  240499.0L;
    constexpr long double  b8 =          170775.0L /                    5593.0L;
    constexpr long double  b6 =         1138500.0L /                   39151.0L;
    constexpr long double  b4 =           14850.0L /                     833.0L;
    constexpr long double  b2 =             108.0L /                      17.0L;

    long double            x2 = a_x * a_x;
    return
      (((((((((((((a27 * x2 + a25) * x2 + a23) * x2 + a21) * x2 + a19)
                       * x2 + a17) * x2 + a15) * x2 + a13) * x2 + a11)
                       * x2 +  a9) * x2 +  a7) * x2 +  a5) * x2 +  a3)
                       * x2 + 1.0L)
      * a_x /
       ((((((((((((b24 * x2 + b22) * x2 + b20) * x2 + b18) * x2 + b16)
                       * x2 + b14) * x2 + b12) * x2 + b10) * x2 + b8)
                       * x2 +  b6) * x2 +  b4) * x2 +  b2) * x2 + 1.0L);
  }

  //=========================================================================//
  // "SqRtPade":                                                             //
  //=========================================================================//
  // Pade approximants for sqrt(x), 1/2 <= x < 1, centered at x = 3/4:
  // Altermatively, we could use approximants centered at x=1, but that would
  // be disadvatageous (similar to "LogPade"):
  //
  template<typename F>
  constexpr F SqRtPade(F a_x);

  //-------------------------------------------------------------------------//
  // SqRtPade<float>: (4,3)th-Order Pade Approximant:                        //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-23) =~ 1.2e-7,
  // the abs precision of the (4,3)th-order approximant is ~1.4e-8:
  //
  template<>
  constexpr float SqRtPade<float>(float a_x)
  {
    assert(0.5f <= a_x && a_x < 1.0f);

    constexpr double A  = SqRt3<double> / 48.0;
    constexpr float  a4 = float(A *   256.0);
    constexpr float  a3 = float(A *  5376.0);
    constexpr float  a2 = float(A * 10080.0);
    constexpr float  a1 = float(A *  3024.0);
    constexpr float  a0 = float(A *    81.0);

    constexpr float  b3 =  64.0f;
    constexpr float  b2 = 336.0f;
    constexpr float  b1 = 252.0f;
    constexpr float  b0 =  27.0f;
    return
      ((((a4 * a_x + a3) * a_x + a2) * a_x + a1) * a_x + a0) /
       (((b3 * a_x + b2) * a_x + b1) * a_x + b0);
  }

  //-------------------------------------------------------------------------//
  // SqRtPade<double>: 8th-Order Pade Approximant:                           //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-52) =~ 2.2e-16,
  // the abs precision of the 8th-order approximation is ~1.5e-17:
  //
  template<>
  constexpr double SqRtPade<double>(double a_x)
  {
    assert(0.5 <= a_x && a_x < 1.0);

    constexpr double A  = SqRt3<double> / 2.0;
    constexpr double a8 = A *   1114112.0;
    constexpr double a7 = A *  33423360.0;
    constexpr double a6 = A * 228114432.0;
    constexpr double a5 = A * 537698304.0;
    constexpr double a4 = A * 504092160.0;
    constexpr double a3 = A * 192471552.0;
    constexpr double a2 = A *  27760320.0;
    constexpr double a1 = A *   1189728.0;
    constexpr double a0 = A *      6561.0;

    constexpr double b8 =     65536.0;
    constexpr double b7 =   6684672.0;
    constexpr double b6 =  87736320.0;
    constexpr double b5 = 342171648.0;
    constexpr double b4 = 504092160.0;
    constexpr double b3 = 302455296.0;
    constexpr double b2 =  72176832.0;
    constexpr double b1 =   5948640.0;
    constexpr double b0 =    111537.0;
    return
      ((((((((a8 * a_x + a7) * a_x + a6) * a_x + a5) * a_x + a4) * a_x + a3)
                 * a_x + a2) * a_x + a1) * a_x + a0) /
      ((((((((b8 * a_x + b7) * a_x + b6) * a_x + b5) * a_x + b4) * a_x + b3)
                 * a_x + b2) * a_x + b1) * a_x + b0);
  }

  //-------------------------------------------------------------------------//
  // SqRtPade<long double>: (10,9)th-Order Pade Approximant:                 //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the (10,9)th-order approximant is ~2.5e-20:
  //
  template<>
  constexpr long double SqRtPade<long double>(long double a_x)
  {
    assert(0.5L <= a_x && a_x < 1.0L);

    constexpr long double A   = SqRt3<long double> / 24.0L;
    constexpr long double a10 = A * 1048576.0L;
    constexpr long double a9  = A * 149422080.0L;
    constexpr long double a8  = A * 2857697280.0L;
    constexpr long double a7  = A * 17146183680.0L;
    constexpr long double a6  = A * 41793822720.0L;
    constexpr long double a5  = A * 45973204992.0L;
    constexpr long double a4  = A * 23509025280.0L;
    constexpr long double a3  = A * 5425159680.0L;
    constexpr long double a2  = A * 508608720.0L;
    constexpr long double a1  = A * 14959080.0L;
    constexpr long double a0  = A * 59049.0L;

    constexpr long double b9  =    1310720.0L;
    constexpr long double b8  =   56033280.0L;
    constexpr long double b7  =  571539456.0L;
    constexpr long double b6  = 2143272960.0L;
    constexpr long double b5  = 3482818560.0L;
    constexpr long double b4  = 2612113920.0L;
    constexpr long double b3  =  904193280.0L;
    constexpr long double b2  =  135628992.0L;
    constexpr long double b1  =    7479540.0L;
    constexpr long double b0  =      98415.0L;
    return
      ((((((((((a10 * a_x + a9) * a_x + a8) * a_x + a7) * a_x + a6) * a_x + a5)
                    * a_x + a4) * a_x + a3) * a_x + a2) * a_x + a1) * a_x + a0)
      / (((((((((b9 * a_x + b8) * a_x + b7) * a_x + b6) * a_x + b5) * a_x + b4)
                    * a_x + b3) * a_x + b2) * a_x + b1) * a_x + b0);
  }

  //=========================================================================//
  // "CbRtPade":                                                             //
  //=========================================================================//
  // Pade approximants for sqrt(x), 1/2 <= x < 1, centered at x = 3/4:
  // Altermatively, we could use approximants centered at x=1, but that would
  // be disadvatageous (similar to "LogPade" and "SqRt"):
  //
  template<typename F>
  constexpr F CbRtPade(F a_x);

  //-------------------------------------------------------------------------//
  // CbRtPade<float>: (4,3)th-Order Pade Approximant:                        //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-23) =~ 1.2e-7,
  // the abs precision of the (4,3)th-order approximant is ~1.4e-8:
  //
  template<>
  constexpr float CbRtPade<float>(float a_x)
  {
    assert(0.5f <= a_x && a_x < 1.0f);

    constexpr double A  = CbRt48<double> / 72.0;
    constexpr float  a4 = float(A *   896.0);
    constexpr float  a3 = float(A * 29568.0);
    constexpr float  a2 = float(A * 66528.0);
    constexpr float  a1 = float(A * 23760.0);
    constexpr float  a0 = float(A *   891.0);

    constexpr float  b3 =  704.0f;
    constexpr float  b2 = 3168.0f;
    constexpr float  b1 = 2079.0f;
    constexpr float  b0 =  189.0f;
    return
      ((((a4 * a_x + a3) * a_x + a2) * a_x + a1) * a_x + a0) /
       (((b3 * a_x + b2) * a_x + b1) * a_x + b0);
  }

  //-------------------------------------------------------------------------//
  // CbRtPade<double>: 8th-Order Pade Approximant:                           //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-52) =~ 2.2e-16,
  // the abs precision of the 8th-order approximation is ~1.5e-17:
  //
  template<>
  constexpr double CbRtPade<double>(double a_x)
  {
    assert(0.5 <= a_x && a_x < 1.0);

    constexpr double A  = CbRt48<double> / 8.0;
    constexpr double a8 = A *    161873920.0;
    constexpr double a7 = A *   5584650240.0;
    constexpr double a6 = A *  41884876800.0;
    constexpr double a5 = A * 106806435840.0;
    constexpr double a4 = A * 107833420800.0;
    constexpr double a3 = A *  44481286080.0;
    constexpr double a2 = A *   7023360960.0;
    constexpr double a1 = A *    342046800.0;
    constexpr double a0 = A *      2565351.0;

    constexpr double b8 =    12812288.0;
    constexpr double b7 =   960921600.0;
    constexpr double b6 = 11098644480.0;
    constexpr double b5 = 39538920960.0;
    constexpr double b4 = 53916710400.0;
    constexpr double b3 = 30039310080.0;
    constexpr double b2 =  6626318400.0;
    constexpr double b1 =   496973880.0;
    constexpr double b0 =     8102835.0;
    return
      ((((((((a8 * a_x + a7) * a_x + a6) * a_x + a5) * a_x + a4) * a_x + a3)
                 * a_x + a2) * a_x + a1) * a_x + a0) /
      ((((((((b8 * a_x + b7) * a_x + b6) * a_x + b5) * a_x + b4) * a_x + b3)
                 * a_x + b2) * a_x + b1) * a_x + b0);
  }

  //-------------------------------------------------------------------------//
  // CbRtPade<long double>: (10,9)th-Order Pade Approximant:                 //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-63) =~ 1.1e-19,
  // the abs precision of the (10,9)th-order approximant is ~1.5e-20:
  //
  template<>
  constexpr long double CbRtPade<long double>(long double a_x)
  {
    assert(0.5L <= a_x && a_x < 1.0L);

    constexpr long double A   = CbRt48<long double> / 288.0L;
    constexpr long double a10 = A * 2789212160.0L;
    constexpr long double a9  = A * 606653644800.0L;
    constexpr long double a8  = A * 13308464332800.0L;
    constexpr long double a7  = A * 87455622758400.0L;
    constexpr long double a6  = A * 229571009740800.0L;
    constexpr long double a5  = A * 270187419156480.0L;
    constexpr long double a4  = A * 147758744851200.0L;
    constexpr long double a3  = A * 36661944211200.0L;
    constexpr long double a2  = A * 3749517021600.0L;
    constexpr long double a1  = A * 124983900720.0L;
    constexpr long double a0  = A * 669556611.0L;

    constexpr long double b9  =    1857781760.0L;
    constexpr long double b8  =   70224150528.0L;
    constexpr long double b7  =  658351411200.0L;
    constexpr long double b6  = 2304229939200.0L;
    constexpr long double b5  = 3518065353600.0L;
    constexpr long double b4  = 2483340249600.0L;
    constexpr long double b3  =  807085581120.0L;
    constexpr long double b2  =  112791463200.0L;
    constexpr long double b1  =    5693799825.0L;
    constexpr long double b0  =      65445975.0L;
    return
      ((((((((((a10 * a_x + a9) * a_x + a8) * a_x + a7) * a_x + a6) * a_x + a5)
                    * a_x + a4) * a_x + a3) * a_x + a2) * a_x + a1) * a_x + a0)
      / (((((((((b9 * a_x + b8) * a_x + b7) * a_x + b6) * a_x + b5) * a_x + b4)
                    * a_x + b3) * a_x + b2) * a_x + b1) * a_x + b0);
  }

  //=========================================================================//
  // Elementary Mathematical Functions for any real types:                   //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "Exp" for an arbitrary real arg:                                        //
  //-------------------------------------------------------------------------//
  // XXX: In CLang <= 19,  this function  is STILL NOT "constexpr", because
  // "std::modf" and "std::ldexp" are not -- so using it in the "constexpr"
  // context will fail:
  //
  template<typename F>
  constexpr F Exp(F a_x)
  {
    // Complex "Exp" has a separate specific implementation:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    // Special Cases:
    if (std::isnan(a_x))
      return NaN<F>;
    if (std::isinf(a_x))
      return (a_x > 0) ? Inf<F> : F(0);

    // Generic Case:
    // First, change the base from E to 2, to make it easier to represent the
    // result in the IEEE 754 format: exp(x) = 2^y, where y = x * Log2E:
    F y = a_x * Log2E<F>;

    // Get the integral and fractional part of "y":   y = intgY + fracY:
    F intgY = NaN<F>;
    F fracY = std::modf (y, &intgY);
    assert(std::isfinite(intgY) && Abs(fracY) < F(1.0));

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
    assert(Abs(fracY) <= F(0.5));

    // For the fractional part, get back to Base-E; it can only become smaller
    // in abs value as a result:
    F f = fracY * Ln2<F>;
    assert(Abs(f) <  F(0.5));

    // Use the Pade approximant for exp(f) around f=0:
    F res = ExpPade<F>(f);

    // Multiply the result by 2^n (can still get 0 or +oo at this moment):
    return std::ldexp(res, n);
# else
    // GCC, and NOT forcing the Pade approximants method:
    return std::exp(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "Log" (Natural Logarithm) for an arbitrary real arg:                    //
  //-------------------------------------------------------------------------//
  // XXX: In CLang <= 19,  this function  is STILL NOT "constexpr", because
  // "std::modf" and "std::ldexp" are not -- so using it in the "constexpr"
  // context will fail:
  //
  template<typename F>
  constexpr F Log(F a_x)
  {
    // Complex "Log" has a separate specific implementation:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    // Special Cases:
    if (std::isnan(a_x) || a_x < 0)
      return NaN<F>;

    if (a_x == 0)
      return -Inf<F>;

    if (a_x == F(1.0))
      return F(0.0);

    if (std::isinf(a_x))
    {
      assert(a_x > 0);
      return  Inf<F>;
    }
    // Generic Case:
    // Get the Base-2 exponent and the normalised fractional part of "a_x":
    assert(std::isfinite(a_x) && a_x > 0);
    int e2X   = INT_MIN;
    F   fracX = std::frexp(a_x, &e2X);
    assert(0.5 <= fracX && fracX < F(1) && e2X != INT_MIN);

    // Then log(x) = log(2^e2X  * fracX) = e2X * log(2) + log(fracX),
    // so expand     log(fracX) in the vicinity of fracX = 3/4:
    F logFX = LogPade<F>(fracX);
    assert(logFX < 0);

    // From "e2X" (which is actually the main part of the Base-2 log of
    // "a_x"), get back to the natural log:
    return F(e2X) * Ln2<F> + logFX;
# else
    // GCC, and NOT forcing the Pade approximants method:
    return std::log(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // Generic "Pow" for an arbitrary real arg:                                //
  //-------------------------------------------------------------------------//
  // XXX: In CLang <= 19,  this function  is STILL NOT "constexpr",  because
  // "Exp" is not either -- so using it in the "constexpr" context will fail:
  //
  template<typename F>
  constexpr F Pow(F a_x, F a_y)
  {
    // Complex "Pow" is not implemented yet. Also, for this generic implementa-
    // tion,   "a_x" is required to be positive:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    assert(a_x > 0);
    return Exp<F>(a_y * Log<F>(a_x));
  }

  //-------------------------------------------------------------------------//
  // "Cos" for an arbitrary real arg:                                        //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F Cos(F a_x)
  {
    // Complex "Cos" has a separate specific implementation:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    if (!std::isfinite(a_x))
      return NaN<F>;

    // First, normalise the arg to the interval [0; +00):
    a_x = Abs(a_x);

    // Then normalise it to the interval [0..2*Pi) (full period):
    a_x = FMod(a_x, TwoPi<F>);

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
    assert(Abs(res) < F(1.0) + Eps<F>);
    return res;
# else
    // GCC, and NOT forcing the Pade approximants method:
    return std::cos(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "Sin" for an arbitrary real arg:                                        //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F Sin(F a_x)
  {
    // Complex "Sin" has a separate specific implementation:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    if (!std::isfinite(a_x))
      return NaN<F>;

    // First, normalise the arg to the interval [0; +00):
    bool  chSgn =  (a_x < 0);
    if (chSgn)
      a_x = - a_x;

    // Then normalise it to the interval [0..2*Pi) (full period):
    a_x = FMod(a_x, TwoPi<F>);

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
    assert(Abs(res) < F(1.0) + Eps<F>);
    return res;
# else
    // GCC, and NOT forcing the Pade approximants method:
    return std::sin(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "Tan" for an arbitrary real arg:                                        //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F Tan(F a_x)
  {
    // Complex "Tan" has a separate specific implementation:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    if (!std::isfinite(a_x))
      return NaN<F>;

    // First, normalise the arg to the interval [0; +00):
    bool  chSgn =  (a_x < 0);
    if (chSgn)
      a_x = - a_x;

    // Then normalise it to the interval [0..Pi) (full period):
    a_x   = FMod(a_x, Pi<F>);

    // Then to the interval [0..Pi/2):
    if (a_x > Pi_2<F>)
    {
      a_x   = Pi<F> - a_x;
      chSgn = !chSgn;
    }

    // Finally, to the interval [0..Pi/4], and compute the function:
    F res = NaN<F>;
    if (a_x <= Pi_4<F>)
      // XXX: For the moment, we do not have a separate Pade Approximant for
      // Tan, so use the ratio. This may result in a larger error, but anyway,
      // the error may be unbounded as x -> Pi/2:
      res = SinPade<F>(a_x) / CosPade<F>(a_x);
    else
    {
      F y = Pi_2<F> - a_x;
      // Tan(x) = CoTan(y):
      res = CosPade<F>(y)   / SinPade<F>(y);
    }
    return chSgn ? (-res) : res;
# else
    // Use the standard impl:
    return std::tan(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "ATan" for an arbitrary real arg:                                       //
  //-------------------------------------------------------------------------//
  // NB: This function is ALWAYS "constexpr", even in CLang:
  //
  template<typename F>
  constexpr F ATan (F a_x)
  { 
    // Complex "ATan" is not implemented yet:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
  
# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    if (a_x == F(0.0))
      return   F(0.0);

    bool chSgn = false;
    if (a_x <  F(0.0))
    {
      a_x = - a_x;
      chSgn = true;
    }
    if (std::isinf(a_x))
      return chSgn ? (-Pi_2<F>) : Pi_2<F>;

    // Then normalise the arg to the interval (0..1]:
    bool inv = false;
    if (a_x > F(1.0))
    {
      a_x = F(1.0) / a_x;
      inv = true;
    }
    // Now use the Pade Approximant:
    F   res = ATanPade<F>(a_x);

    if (inv)
      res = Pi_2<F> - res;
    if (chSgn)
      res = - res;
    return res;
# else
    // Use the standard impl:
    return std::atan(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "ATan2" for an arbitrary real arg:                                      //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F ATan2(F a_y, F a_x)
  {
    // The args must be real:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    // "a_y" has the sign of Sin, "a_x" has the sign of Cos:
    // XXX: this impl is NOT completely conformant to that of std::atan2; some
    // edge cases are handled differently here:
    return
      (a_x >  F(0.0))                     // 1st and 4th Quadrants
      ? ATan<F>(a_y  /  a_x)
      :
      (a_x <  F(0.0))
      ? ((a_y >= F(0.0))
         ? ATan<F>(a_y  / a_x) + Pi<F>    // 2nd Quadrant
         : ATan<F>(a_y  / a_x) - Pi<F>)   // 3rd Quadrant
      :
      ((a_y > F(0.0))
       ?  Pi_2<F>                         // Pos Y Axis
       :
       (a_y < F(0.0))                     // Neg Y Axis
       ? -Pi_2<F>
       : NaN<F>);                         // 0/0 is UNDEFINED
# else
    // Use the standard impl:
    return std::atan2(a_y, a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "SqRt" (Square Root) for an arbitrary real arg:                         //
  //-------------------------------------------------------------------------//
  // NB: This function is ALWAYS "constexpr", even in CLang:
  //
  template<typename F>
  constexpr F SqRt(F a_x)
  {
    // Complex "SqRt" has a separate specific implementation:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    // Special Cases:
    if (!(a_x >= 0))
      return NaN<F>;

    if (a_x == 0 || a_x == F(1.0))
      return a_x;

    if (std::isinf(a_x))
    {
      assert(a_x > 0);
      return  Inf<F>;
    }
    // Generic Case:
    //
#   if (!defined(__clang__))
    // Get the Base-2 exponent and the normalised fractional part of "a_x":
    assert(std::isfinite(a_x) && a_x > 0);
    int e2X   = INT_MIN;
    F   fracX = std::frexp(a_x, &e2X);
    assert(0.5 <= fracX && fracX < F(1) && e2X != INT_MIN);

    // Get the Pade approximant of SqRt(fracX):
    F sqrtFX = SqRtPade<F>(fracX);
    assert(0 < sqrtFX && sqrtFX < F(1.0) + Eps<F>);

    // Then the exponent of the result is e2X/2:
    int n = e2X / 2;
    if (e2X % 2 != 0)
    {
      if (e2X > 0)
        sqrtFX *= SqRt2<F>;
      else
        sqrtFX /= SqRt2<F>;
    }
    return std::ldexp(sqrtFX, n);
#   else
    // XXX: In CLang <= 19, "std::frexp" is not "constexpr", so another method
    // is required: Use Hayley iterations which always converge:
    // Initial approximation of the root from above: Using the ineqiality bet-
    // ween the Arithmetic and Geometric Means (the other operand is 1):
    //
    F y = F(0.5) * (a_x + F(1.0));

    while (true)
    {
      // The following normalisation tries to avoid Infinities and NaNs:
      F x_y2 = (a_x / y) / y;
      F c    = (F(1.0) + F(3.0) * x_y2) / (F(3.0) + x_y2);
      y *= c;

      // "ry" is a relative correction applied to "y":
      F ry = c - F(1.0);
      if (ry < 0)
        ry = - ry;
      if (ry < F(50.0) * Eps<F>)
        break;
    }
    return y;
#   endif
# else
    // GCC, and NOT forcing the Pade approximants method:
    return std::sqrt(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "CbRt" for an arbitrary real arg:                                       //
  //-------------------------------------------------------------------------//
  // NB: This function is ALWAYS "constexpr", even in CLang:
  //
  template<typename F>
  constexpr F CbRt (F a_x)
  {
    // Complex "SqRt" has a separate specific implementation:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    if (std::isinf(a_x))
      return (a_x > 0) ? Inf<F> : -Inf<F>;
    else
    if (a_x == 0)
      return F(0.0);
    // NaN is not considered specifically...

    // Generic Case:
    bool chSgn = (a_x < 0);
    if (chSgn)
      a_x = - a_x;
    assert(a_x > 0);

#   if (!defined(__clang__))
    // Get the Base-2 exponent and the normalised fractional part of "a_x":
    assert(std::isfinite(a_x) && a_x > 0);
    int e2X   = INT_MIN;
    F   fracX = std::frexp(a_x, &e2X);
    assert(0.5 <= fracX && fracX < F(1) && e2X != INT_MIN);

    // Get the Pade approximant of CbRt(fracX):
    F cbrtFX = CbRtPade<F>(fracX);
    assert(0 < cbrtFX && cbrtFX < F(1.0) + Eps<F>);

    // Then the exponent of the result is e2X/3:
    int n = e2X / 3;
    int r = e2X % 3;  // 0, +-1, +-2
    switch (r)
    {
      case 1:
        assert(e2X > 0 && n >= 0);
        cbrtFX *= CbRt2<F>;
        break;
      case -1:
        assert(e2X < 0 && n <= 0);
        cbrtFX /= CbRt2<F>;
        break;
      case 2:
        assert(e2X > 0 && n >= 0);
        cbrtFX *= CbRt4<F>;
        break;
      case -2:
        assert(e2X < 0 && n <= 0);
        cbrtFX /= CbRt4<F>;
        break;
      case 0 :
        break;          // No adjustments
      default:
        assert(false);  // Cannot happen
    }
    F y = std::ldexp(cbrtFX, n);
#   else
    // XXX: In CLang <= 19, "std::frexp" is not "constexpr", so another method
    // is required: Use Hayley iterations which always converge:
    // Initial approximation of the root from above: Using the ineqiality bet-
    // ween the Arithmetic and Geometric Means (the other operands are 1, 1):
    //
    F y = (a_x + F(2.0)) /  F(3.0);

    while (true)
    {
      // The following normalisation tries to avoid Infinities and NaNs:
      F x_y3 = ((a_x / y) / y) / y;
      F c    = (F(1.0) + F(2.0) * x_y3) / (F(2.0) + x_y3);
      y *= c;

      // "ry" is a relative correction applied to "y":
      F ry = c - F(1.0);
      if (ry < 0)
        ry = - ry;
      if (ry < F(50.0) * Eps<F>)
        break;
    }
#   endif
    if (chSgn)
      y = -y;
    return y;
# else
    // GCC, and NOT forcing the Pade approximants method:
    return std::cbrt (a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "ASin" for an arbitrary real arg:                                       //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F ASin (F a_x)
  {
    // Complex "ASin" is not implemented yet:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    return
      (-F(1.0) < a_x && a_x < F(1.0))
      ? ATan(a_x / SqRt(F(1.0) - a_x * a_x))
      :
      (a_x == -F(1.0))
      ? -Pi_2<F>
      :
      (a_x ==  F(1.0))
      ?  Pi_2<F>
      :  NaN<F>;
# else
    // Use the standard impl:
    return std::asin(a_x);
# endif
  }

  //-------------------------------------------------------------------------//
  // "ACos" for an arbitrary real arg:                                       //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F ACos (F a_x)
  {
    // Complex "ACos" is not implemented yet:
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);

# if (defined(__clang__) || DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL)
    return
      (F(0.0) < a_x   && a_x <= F(1.0))
      ? ATan(SqRt(F(1.0) - a_x * a_x) / a_x)
      :
      (a_x == F(0.0))
      ? Pi_2<F>
      :
      (-F(1.0) <= a_x && a_x <  F(0.0))
      ? ATan(SqRt(F(1.0) - a_x * a_x) / a_x) + Pi<F>
      : NaN<F>;
# else
    // Use the standard impl:
    return std::acos(a_x);
# endif
  }

  //=========================================================================//
  // Hyperbolic and Inverse Hyperbolic Functions:                            //
  //=========================================================================//
  template<typename F>
  constexpr F SinH (F a_x)
  {
    static_assert(std::is_floating_point_v<F>);
    F ex = Exp(a_x);
    return (ex - F(1.0) / ex) / F(2.0);
  }

  template<typename F>
  constexpr F CosH (F a_x)
  {
    static_assert(std::is_floating_point_v<F>);
    F ex = Exp(a_x);
    return (ex + F(1.0) / ex) / F(2.0);
  }

  template<typename F>
  constexpr F TanH (F a_x)
  {
    static_assert(std::is_floating_point_v<F>);
    F ex2 = Exp(2.0 * a_x);
    return (ex2 - F(1.0)) / (ex2 + F(1.0));
  }

  // XXX: The Inverse Hyperbolic Functions are currently implemented for Real
  // args only:
  //
  template<typename F>
  constexpr F ASinH(F a_x)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    return Log(a_x + SqRt(Sqr(a_x) + F(1.0)));
  }

  template<typename F>
  constexpr F ACosH(F a_x)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    assert(a_x >= F(1.0));
    // Return the non-negative root:
    return Log(a_x + SqRt(Sqr(a_x) - F(1.0)));
  }

  template<typename F>
  constexpr F ATanH(F a_x)
  {
    static_assert(std::is_floating_point_v<F> && !IsComplex<F>);
    assert(Abs(a_x) < F(1.0));
    return F(0.5) * Log((F(1.0) + a_x) / (F(1.0) - a_x));
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
  // XXX: Not really "constexpr" in CLang <= 19:
  //
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
  // XXX: Not really "constexpr" in CLang <= 19:
  //
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
  // XXX: Not really "constexpr" in CLang <= 19:
  //
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
  // Returns (cos(z), sin(z)). XXX: Not really "constexpr" in CLang <= 19:
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
  // Complex "Tan":                                                          //
  //-------------------------------------------------------------------------//
  // NB: This function is ALWAYS "constexpr", even in CLang:
  //
  template<typename T>
  constexpr std::complex<T> Tan(std::complex<T> a_z)
  {
    // Complex arg: Compute Cos and Sin together:
    auto   cs = CosSin(a_z);
    return std::get<1>(cs) / std::get<0>(cs);
  }

  //-------------------------------------------------------------------------//
  // Complex "Pow":                                                          //
  //-------------------------------------------------------------------------//
  // XXX: Not really "constexpr" in CLang <= 19:
  //
  template<typename T>
  constexpr  std::complex<T> Pow(std::complex<T> a_z, T a_p)
    { return std::pow(a_z, a_p); }

  //-------------------------------------------------------------------------//
  // Complex "SqRt":                                                         //
  //-------------------------------------------------------------------------//
  // XXX: Also a placeholder as yet, NOT "constexpr" in CLang <= 19:
  //
  template<typename T>
  constexpr  std::complex<T> SqRt(std::complex<T> a_z)
    { return std::sqrt(a_z); }

  //-------------------------------------------------------------------------//
  // Complex "CbRt":                                                         //
  //-------------------------------------------------------------------------//
  // XXX: Also a placeholder as yet, NOT "constexpr" in CLang <= 19:
  //
  template<typename T>
  constexpr  std::complex<T> CbRt(std::complex<T> a_z)
    { return std::pow(a_z, T(1)/T(3)); }
}
