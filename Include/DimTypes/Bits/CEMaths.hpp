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
  // Eps<float> = 2^(-24) =~ 6e-8,
  // the abs precision of the 4th-order approximant is ~1e-10 (however, the
  // 3rd-order one would be slightly insuffient):
  //
  template<>
  constexpr float ExpPade<float>(float a_x)
  {
    assert(std::fabs(a_x) <= 0.5f);

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
    return (se + so) / (se - so);
  }

  //-------------------------------------------------------------------------//
  // ExpPade<double>: 6th-Order Pade Approximant:                            //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-53) =~ 1e-16,
  // the abs precision of the 6th-order approximation is ~4e-17:
  //
  template<>
  constexpr double ExpPade<double>(double a_x)
  {
    assert(std::fabs(a_x) <= 0.5);

    constexpr double a1 = 0.5;
    constexpr double a2 = 5.0 /     44.0;
    constexpr double a3 = 1.0 /     66.0;
    constexpr double a4 = 1.0 /    792.0;
    constexpr double a5 = 1.0 /  15840.0;
    constexpr double a6 = 1.0 / 665280.0;

    double x2  = a_x * a_x;
    double x3  = x2  * a_x;
    double x4  = x2  * x2;
    double x5  = x4  * a_x;
    double x6  = x3  * x3;

    double t1  = a1  * a_x;
    double t2  = a2  * x2;
    double t3  = a3  * x3;
    double t4  = a4  * x4;
    double t5  = a5  * x5;
    double t6  = a6  * x6;

    double se  = t6  + t4  + t2 + 1.0;
    double so  = t5  + t3  + t1;
    return (se + so) / (se - so);
  }

  //-------------------------------------------------------------------------//
  // ExpPade<long double>: 7th-Order Pade Approximant:                       //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-64) =~ 5e-20,
  // the abs precision of the 7th-order approximant is ~1e-20:
  //
  template<>
  constexpr long double ExpPade<long double>(long double a_x)
  {
    assert(std::fabs(a_x) <= 0.5L);

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
    return (se + so) / (se - so);
  }

  //=========================================================================//
  // "LogPade":                                                              //
  //=========================================================================//
  // Pade approximants for log(x), 1/2 <= x < 1, centered at x = 3/4:
  //
  template<typename F>
  F LogPade(F a_x);

  //-------------------------------------------------------------------------//
  // LogPade<float>: 4th-Order Pade Approximant:                             //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-24) =~ 6e-8,
  // the abs precision of the 4th-order approximant is ~6e-9:
  //
  template<>
  constexpr float LogPade<float>(float a_x)
  {
    assert(0.5f <= a_x && a_x < 1.0f);
    float y = a_x - 0.75f;

    constexpr float  b4 =   384.0f;
    constexpr float  b3 =  5760.0f;
    constexpr float  b2 = 19440.0f;
    constexpr float  b1 = 22680.0f;
    constexpr float  b0 =  8505.0f;

    constexpr double A  = -0.28768207245178093;         // log(3/4)
    constexpr float  a4 = float(  384.0 * A +  1600.0);
    constexpr float  a3 = float( 5760.0 * A + 12480.0);
    constexpr float  a2 = float(19440.0 * A + 22680.0);
    constexpr float  a1 = float(22680.0 * A + 11340.0);
    constexpr float  a0 = float( 8505.0 * A);
    return
      ((((a4 * y + a3) * y + a2) * y + a1) + a0) /
      ((((b4 * y + b3) * y + b2) * y + b1) + b0);
  }

  //-------------------------------------------------------------------------//
  // LogPade<double>: 8th-Order Pade Approximant:                            //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-53) =~ 1e-16,
  // the abs precision of the 8th-order approximation is ~7e-17:
  //
  template<>
  constexpr double LogPade<double>(double a_x)
  {
    assert(0.5 <= a_x && a_x < 1.0);
    double y = a_x - 0.75;

    constexpr double b8 =     1146880.0;
    constexpr double b7 =    61931520.0;
    constexpr double b6 =   812851200.0;
    constexpr double b5 =  4470681600.0;
    constexpr double b4 = 12573792000.0;
    constexpr double b3 = 19615115520.0;
    constexpr double b2 = 17163226080.0;
    constexpr double b1 =  7881073200.0;
    constexpr double b0 =  1477701225.0;

    constexpr double A  =  -0.28768207245178093;          // log(3/4)
    constexpr double a8 =  b8 * A +     6234112.0;
    constexpr double a7 =  b7 * A +   212779008.0;
    constexpr double a6 =  b6 * A +  1979873280.0;
    constexpr double a5 =  b5 * A +  7908848640.0;
    constexpr double a4 =  b4 * A + 15956740800.0;
    constexpr double a3 =  b3 * A + 17046469440.0;
    constexpr double a2 =  b2 * A +  9194585400.0;
    constexpr double a1 =  b1 * A +  1970268300.0;
    constexpr double a0 =  b0 * A;
    return
      ((((((((a8 * y + a7) * y + a6) * y + a5) * y + a4) * y + a3) * y + a2)
                 * y + a1) * y + a0) /
      ((((((((b8 * y + b7) * y + b6) * y + b5) * y + b4) * y + b3) * y + a2)
                 * y + b1) * y + b0);
  }

  //-------------------------------------------------------------------------//
  // LogPade<long double>: 10th-Order Pade Approximant:                      //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-64) =~ 5e-20,
  // the abs precision of the 10th-order approximant is ~8e-21:
  //
  template<>
  constexpr long double LogPade<long double>(long double a_x)
  {
    assert(0.5L  <= a_x && a_x < 1.0L);
    long double y = a_x - 0.75L;

    constexpr long double b10 =       82575360.0L;
    constexpr long double b9  =     6812467200.0L;
    constexpr long double b8  =   137952460800.0L;
    constexpr long double b7  =  1195587993600.0L;
    constexpr long double b6  =  5492232345600.0L;
    constexpr long double b5  = 14829027333120.0L;
    constexpr long double b4  = 24715045555200.0L;
    constexpr long double b3  = 25723822924800.0L;
    constexpr long double b2  = 16278356694600.0L;
    constexpr long double b1  =  5727569948100.0L;
    constexpr long double b0  =   859135492215.0L;

    constexpr long double A   = -0.287682072451780927439L;   // log(3/4)
    constexpr long double a10 = b10 * A +      483721216.0L;
    constexpr long double a9  = b9  * A +    26282065920.0L;
    constexpr long double a8  = b8  * A +   394259374080.0L;
    constexpr long double a7  = b7  * A +  2619855912960.0L;
    constexpr long double a6  = b6  * A +  9288846927360.0L;
    constexpr long double a5  = b5  * A + 19148275770624.0L;
    constexpr long double a4  = b4  * A + 23675444432640.0L;
    constexpr long double a3  = b3  * A + 17292125410560.0L;
    constexpr long double a2  = b2  * A +  6873083937720.0L;
    constexpr long double a1  = b1  * A +  1145513989620.0L;
    constexpr long double a0  = b0  * A;
    return
      ((((((((((a10 * y + a9) * y + a8) * y + a7) * y + a6) * y + a5) * y + a4)
                    * y + a3) * y + a2) * y + a1) * y + a0)
      /
      ((((((((((b10 * y + b9) * y + b8) * y + b7) * y + b6) * y + b5) * y + b4)
                    * y + b3) * y + b2) * y + b1) * y + b0);
  }

  //=========================================================================//
  // "CosPade":                                                              //
  //=========================================================================//
  // Assumes 0 <= x <= Pi/4, approximants are cenetred at Pi/8:
  //
  template<typename F>
  F CosPade(F a_x);

  //-------------------------------------------------------------------------//
  // CosPade<float>: 4th-Order Pade approximant:                             //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-24) =~ 6e-8,
  // the abs precision of the 4th-order approximant is ~4e-10 (yet, lower-order
  // approximants are insufficient):
  //
  template<>
  constexpr float CosPade<float>(float a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<float> + 10.0f * Eps<float>);
    float  y  = a_x - Pi_8<float>;

    constexpr float C    = CosPi8<float>;
    constexpr float C2   = C  * C;
    constexpr float C4   = C2 * C2;
    constexpr float S    = SinPi8<float>;
    constexpr float S2   = S  * S;
    constexpr float S4   = S2 * S2;
    constexpr float CS   = C  * S;
    constexpr float C2S2 = CS * CS;

    constexpr float a4 =         C *(32865.0f*C4 + 63720.0f*C2S2 + 30856.0f*S4);
    constexpr float a3 =   20.0f*S * (6573.0f*C4 + 12028.0f*C2S2 +  5456.0f*S4);
    constexpr float a2 = -180.0f*C * (4025.0f*C4 +  7312.0f*C2S2 +  3288.0f*S4);
    constexpr float a1 = -840.0f*S * (1725.0f*C4 +  2956.0f*C2S2 +  1232.0f*S4);
    constexpr float a0 = 1680.0f*C * ( 945.0f*C4 +  1560.0f*C2S2 +   616.0f*S4);

    constexpr float b4 =    1365.0f*C4 +    3300.0f*C2S2 +    1936.0f*S4;
    constexpr float b3 =      20.0f*CS *    (273.0f*C2   +     274.0f*S2);
    constexpr float b2 =   69300.0f*C4 +  132840.0f*C2S2 +   63360.0f*S4;
    constexpr float b1 =     840.0f*CS *    (165.0f*C2   +     164.0f*S2);
    constexpr float b0 = 1587600.0f*C4 + 2620800.0f*C2S2 + 1034880.0f*S4;
    return
      ((((a4 * y + a3) * y + a2) * y + a1) * y + a0) /
      ((((b4 * y + b3) * y + b2) * y + b1) * y + b0);
  }

  //-------------------------------------------------------------------------//
  // CosPade<double>: 8th-Order Pade approximant:                            //
  //-------------------------------------------------------------------------//
  // Eps<double> = 2^(-53) =~ 1e-16,
  // the abs precision of the 7th-order approximation is ~2e-18 (however, the
  // 6th-order one would be insufficient):
  //
  template<>
  constexpr double CosPade(double a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<double> + 100.0 * Eps<double>);
    double y  = a_x - Pi_8<double>;

    constexpr double C    = CosPi8<double>;
    constexpr double C2   = C  * C;
    constexpr double C4   = C2 * C2;
    constexpr double C6   = C4 * C2;
    constexpr double C8   = C4 * C4;

    constexpr double S    = SinPi8<double>;
    constexpr double S2   = S  * S;
    constexpr double S4   = S2 * S2;
    constexpr double S6   = S4 * S2;
    constexpr double S8   = S4 * S4;

    constexpr double C6S2 = C6 * S2;
    constexpr double C4S4 = C4 * S4;
    constexpr double C2S6 = C2 * S6;
    constexpr double C4S2 = C4 * S2;
    constexpr double C2S4 = C2 * S4;
    constexpr double CS   = C  * S;

    constexpr double a7 =
      -  17960237039745.0*C8      -  83207769964704.0*C6S2
      - 141933814031520.0*C4S4    - 106085270418688.0*C2S6
      -  29398989312128.0*S8;

    constexpr double a6 =
      56.0*CS *
      (2565748148535.0*C6         + 7730465927976.0*C4S2 +
       7763690345904.0*C2S4       + 2598972566464.0*S6);

    constexpr double a5 =
       1321096787226216.0*C8      +  7172883472132224.0*C6S2 +
      13620536004096000.0*C4S4    + 11006812103261184.0*C2S6 +
       3238062784072704.0*S8;

    constexpr double a4 =
      -5040.0*CS *
      (2621223784179.0*C6         + 7919753370408.0*C4S2 +
       7975844019696.0*C2S4       + 2677314433472.0*S6);

    constexpr double a3 =
      - 22699576770670800.0*C8    - 166809791770997760.0*C6S2
      -366163436979648000.0*C4S4  - 322696181854126080.0*C2S6
      -100642959875082240.0*S8;

    constexpr double a2 =
      1995840.0*CS *
      (136481341815.0*C6          + 413843892840.0*C4S2  +
       418244830896.0*C2S4        + 140882279872.0*S6);

    constexpr double a1 =
        48235826436557760.0*C8    +  826617033671869440.0*C6S2 +
      2214669463598592000.0*C4S4  + 2142438414048092160.0*C2S6 +
       706150157693460480.0*S8;

    constexpr double a0 =
      -17297280.0*CS *
      ( 39040911063.0*C6          + 118905514056.0*C4S2  +
       120688949808.0*C2S4        +  40824346816.0*S6);

    constexpr double b7 =
      C *
      (156069114201.0*C6          + 462824834472.0*C4S2  +
       457443575568.0*C2S4        + 150687855296.0*S6);

    constexpr double b6 =
      -56.0*S *
      (22295587743.0*C6           + 64703433366.0 *C4S2  +
       62520967656.0*C2S4         + 20113122032.0 *S6);

    constexpr double b5 =
      1512.0*C *
      (13517797293.0*C6           + 40955931126.0 *C4S2  +
       41357942008.0*C2S4         + 13919808176.0 *S6);

    constexpr double b4 =
      -5040.0*S *
      ( 40553391879.0*C6          + 119785261860.0*C4S2  +
       117909181272.0*C2S4        +  38677311296.0*S6);

    constexpr double b3 =
      55440.0*C *
      (25583269257.0 *C6          +  78971217948.0*C4S2  +
       81192490200.0 *C2S4        +  27804541504.0*S6);

    constexpr double b2 =
      -1995840.0*S *
      ( 8527756419.0 *C6          + 25597906866.0 *C4S2  +
       25612284624.0 *C2S4        +  8542134176.0 *S6);

    constexpr double b1 =
      8648640.0*C *
      (5577273009.0  *C6          + 17495867970.0 *C4S2  +
       18260384688.0 *C2S4        +  6341789728.0 *S6);

    constexpr double b0 =
      -17297280.0*S *
      ( 39040911063.0*C6          + 118905514056.0*C4S2  +
       120688949808.0*C2S4        +  40824346816.0*S6);

    return
      (((((((a7 * y + a6) * y + a5) * y + a4) * y + a3) * y + a2) * y + a1)
                * y + a0)   /
      (((((((b7 * y + b6) * y + b5) * y + b4) * y + b3) * y + b2) * y + b1)
                * y + b0);
  }

  //-------------------------------------------------------------------------//
  // CosPade<long double>: 8th-Order Pade approximant:                       //
  //-------------------------------------------------------------------------//
  // Eps<long double> = 2^(-64) =~ 5e-20,
  // the abs precision of the 8th-order approximant is ~1e-22:
  //
  template<>
  constexpr long double CosPade(long double a_x)
  {
    assert(0 <= a_x && a_x < Pi_4<long double> + 100.0L * Eps<long double>);
    long double y   =  a_x - Pi_8<long double>;

    constexpr long double C    = CosPi8<long double>;
    constexpr long double C2   = C  * C;
    constexpr long double C4   = C2 * C2;
    constexpr long double C6   = C4 * C2;
    constexpr long double C8   = C4 * C4;

    constexpr long double S    = SinPi8<long double>;
    constexpr long double S2   = S  * S;
    constexpr long double S4   = S2 * S2;
    constexpr long double S6   = S4 * S2;
    constexpr long double S8   = S4 * S4;

    constexpr long double C6S2 = C6 * S2;
    constexpr long double C4S4 = C4 * S4;
    constexpr long double C2S6 = C2 * S6;
    constexpr long double C4S2 = C4 * S2;
    constexpr long double C2S4 = C2 * S4;

    constexpr long double a8 =
      C*
      ( 297652211504544897.0L*C8        + 1174264331449521312.0L*C6S2 +
       1736890268152061088.0L*C4S4      + 1141596387920467200.0L*C2S6 +
        281318239713382528.0L*S8);

    constexpr long double a7 =
      72.0L*S*
      ( 33072467944949433.0L*C8         + 127323946914556392.0L*C6S2  +
       183542490714637296.0L*C4S4       + 117403012425542080.0L*C2S6  +
        28112000680511744.0L*S8);

    constexpr long double a6 =
      -2520.0L*C*
      ( 20267255425871025.0L*C8         + 78716525642582256.0L*C6S2   +
       114549299020529472.0L*C4S4       + 74018042785709184.0L*C2S6   +
        17918013981890944.0L*S8);

    constexpr long double a5 =
      -55440.0L*S*
      ( 5527433297964825.0L*C8          + 20979200214016200.0L*C6S2   +
       29775043377698160.0L*C4S4        + 18722219280440512.0L*C2S6   +
        4398942818793728.0L*S8);

    constexpr long double a4 =
      166320.0L*C*
      (16130697800718501.0L*C8          + 61507376437867776.0L*C6S2   +
       87744575600843232.0L*C4S4        + 55489812991102464.0L*C2S6   +
       13121916027408512.0L*S8);

    constexpr long double a3 =
      8648640.0L*S*
      (1240822907747577.0L*C8           + 4636795860808488.0L*C6S2    +
       6466327458527472.0L*C4S4         + 3985558949424064.0L*C2S6    +
        915204443957504.0L*S8);

    constexpr long double a2 =
      -8648640.0L*C*
      ( 4806222105259575.0L*C8          + 17941474161593136.0L*C6S2   +
       24991189867060992.0L*C4S4        + 15382845578391168.0L*C2S6   +
        3526907767663744.0L*S8);

    constexpr long double a1 =
      -259459200.0L*S*
      ( 320414807017305.0L*C8           + 1177883440344072.0L*C6S2    +
       1611553286068848.0L*C4S4         +  971115468558016.0L*C2S6    +
        217030815815936.0L*S8);

    constexpr long double a0 =
      518918400.0L*C*
      (167629288667841.0L*C8            +  611142073397856.0L*C6S2    +
       827911697193504.0L*C4S4          +  492914320371456.0L*C2S6    +
       108515407907968.0L*S8);

    constexpr long double b8 =
       542578576637097.0L*C8            + 2427989179575960.0L*C6S2    +
      4029330951462960.0L*C4S4          + 2945008687202112.0L*C2S6    +
       801088338678016.0L*S8;

    constexpr long double b7 =
      72.0L*C*S*
      ( 60286508515233.0L*C6            + 181255851122346.0L*C4S2     +
       181652188684824.0L*C2S4          +  60682846077712.0L*S6);

    constexpr long double b6 =
      144992827721621880.0L*C8          + 601083367536571920.0L*C6S2  +
      932967718803489600.0L*C4S4        + 642656625645525120.0L*C2S6  +
      165779446656983040.0L*S8;

    constexpr long double b5 =
      55440.0L*C*S*
      (15691864472037.0L*C6             + 47087897200812.0L*C4S2      +
       47100196475016.0L*C2S4           + 15704163746240.0L*S6);

    constexpr long double b4 =
       21282423302370550320.0L*C8       + 83492524053917858880.0L*C6S2     +
      122732020270552049280.0L*C4S4     + 80116162718853381120.0L*C2S6     +
       19594243199849472000.0L*S8;

    constexpr long double b3 =
      8648640.0L*C*S*
      ( 9843130620477.0L*C6             + 29486873498382.0L*C4S2      +
       29444353658160.0L*C2S4           +  9800610780256.0L*S6);

    constexpr long double b2 =
       1925676385894920859200.0L*C8     + 7248435129017981458560.0L*C6S2   +
      10189674067987768273920.0L*C4S4   + 6336748537167824455680.0L*C2S6   +
       1469833212303056240640.0L*S8;

    constexpr long double b1 =
      259459200.0L*C*S*
      (14843770318377.0L*C6             + 44400706451640.0L*C4S2      +
       44270108318160.0L*C2S4           + 14713172184896.0L*S6);

    constexpr long double b0 =
       86985922268654183174400.0L*C8    + 317132866900297998950400.0L*C6S2 +
      429618613248937586073600.0L*C4S4  + 255782310464243353190400.0L*C2S6 +
       56310641846950101811200.0L*S8;

    return
      ((((((((a8 * y + a7) * y + a6) * y + a5) * y + a4) * y + a3) * y + a2)
                 * y + a1) * y + a0)   /
      ((((((((b8 * y + b7) * y + b6) * y + b5) * y + b4) * y + b3) * y + b2)
                 * y + b1) * y + b0);
  }

  //=========================================================================//
  // "SinPade":                                                              //
  //=========================================================================//
  // The arg is assumed to be in [0 .. Pi/4], centered at Pi/8:
  //
  template<typename F>
  F SinPade(F a_x);

  //-------------------------------------------------------------------------//
  // SinPade<float>: (5,4)th-Order Pade approximant:                         //
  //-------------------------------------------------------------------------//
  // Eps<float> = 2^(-24) =~ 6e-8,
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
  // Eps<double> = 2^(-53) =~ 1e-16,
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
  // Eps<long double> = 2^(-64) =~ 5e-20,
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
  template<bool ForceApprox=false, typename F>
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
      assert(std::isfinite(intgY) && std::fabs(fracY) < F(1));

      // Check if "intgY" is too large in absolute value:
      if (intgY > F(INT_MAX))
        return Inf<F>;
      if (intgY < F(INT_MIN))
        return F(0);

      // If OK:
      int n = int(intgY);

      // For the fractional part, get back to Base-E, it can only become smaller
      // in abs value as a result:
      F f = fracY * Ln2<F>;
      assert(std::fabs(f) < F(1));

      // Make it within [-0.5 .. 0.5] for better convergence:
      if (f < F(-0.5))
      {
        f += F(1.0);
        --n;
      }
      else
      if (f > F(0.5))
      {
        f -= F(1.0);
        ++n;
      }
      assert(std::fabs(f) <= F(0.5));

      // Use the Pade approximant for exp(f) around f=0:
      F res = ExpPade<F>(f);

      // Multiply the result by 2^n (can still get 0 or +oo):
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
  template<bool ForceApprox=false, typename F>
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
      // so expand     log(fracX) in the vicinity of fracX = 1/4:
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
  template<bool ForceApprox=false, typename F>
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
  template<bool ForceApprox=false, typename F>
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
      return chSgn ? (- res) : res;
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
  template<bool ForceApprox=false, typename F>
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
      return chSgn ? (- res) : res;
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
