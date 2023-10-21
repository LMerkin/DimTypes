// vim:ts=2:et
//===========================================================================//
//                             "ConstExprMaths.hpp":                         //
//      Standard Mathematical Functions as "constexpr"s, for C++ <= 23       //
//===========================================================================//
#pragma  once
#include <cmath>
#include <cassert>
#include <cfloat>

namespace DimTypes
{
  namespace ConstExprMaths
  {
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // Pi:
    template<typename F> constexpr F Pi    = F(M_PIl);
    // 2*Pi:
    template<typename F> constexpr F TwoPi = F(2.0L * M_PIl);
    // Pi/2:
    template<typename F> constexpr F Pi_2  = F(M_PI_2l);
    // Pi/4:
    template<typename F> constexpr F Pi_4  = F(M_PI_4l);

    //=======================================================================//
    // Pade Approximations for Sin and Cos:                                  //
    //=======================================================================//
    // We want all functions below to be "constexpr", but std::sin and std::cos
    // are not "constexpr" funcs for C++ < 26, so we provide our own, reasonably
    // accurate, implementation, using rational Pade approximations.
    // Absolute error is < 2e-10 which is OK for our applications:
    //
    //-----------------------------------------------------------------------//
    // "SinPade": Arg is assumed to be in [0..Pi/4]:                         //
    //-----------------------------------------------------------------------//
    template<typename F>
    constexpr F SinPade(F a_x)
    {
      assert(0.0 <= a_x && a_x < Pi_4<F> + F(100.0 * DBL_EPSILON));

      constexpr F a5 = F( 12671) / F( 4363920);
      constexpr F a3 = F(- 2363) / F(   18183);
      constexpr F b6 = F(   121) / F(16662240);
      constexpr F b4 = F(   601) / F(  872784);
      constexpr F b2 = F(   445) / F(   12122);
      F           x2 = a_x * a_x;
      return a_x * ( (a5 * x2 + a3) * x2 + F(1))
                 / (((b6 * x2 + b4) * x2 + b2) * x2 + F(1));
    }

    //-----------------------------------------------------------------------//
    // "CosPade": Arg is assumed to be in [0..Pi/4]:                         //
    //-----------------------------------------------------------------------//
    template<typename F>
    constexpr F CosPade(F a_x)
    {
      assert(0 <= a_x && a_x < Pi_4<F> + F(100.0 * DBL_EPSILON));
      constexpr F a4 = F( 30257) / F(1577520);
      constexpr F a2 = F(-  425) / F(    939);
      constexpr F b6 = F(    59) / F(3155040);
      constexpr F b4 = F(  1907) / F(1577520);
      constexpr F b2 = F(    89) / F(   1878);
      F           x2 = a_x * a_x;
      return ( (a4 * x2 + a2) * x2 + F(1))  /
             (((b6 * x2 + b4) * x2 + b2) * x2 + F(1));
    }

    //-----------------------------------------------------------------------//
    // "Sin" for an arbitrary arg:                                           //
    //-----------------------------------------------------------------------//
    template<typename F>
    constexpr F Sin(F a_x)
    {
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
      F res = (a_x <= Pi_4<F>) ? SinPade(a_x) : CosPade(Pi_2<F> - a_x);

      // Don't forget the sign:
      return chSgn ? (- res) : res;
    }

    //-----------------------------------------------------------------------//
    // "Cos" for an arbitrary arg:                                           //
    //-----------------------------------------------------------------------//
    template<typename F>
    constexpr F Cos(F a_x)
    {
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
      F res = (a_x <= Pi_4<F>) ? CosPade(a_x) : SinPade(Pi_2<F> - a_x);

      // Don't forget the sign:
      return chSgn ? (- res) : res;
    }
  }

  //=========================================================================//
  // Squre and Cubic Roots by Halley Iterations:                             //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "SqRt":                                                                 //
  //-------------------------------------------------------------------------//
  template<typename F>
  constexpr F SqRt (F a_x)
  {
    if (a_x < 0)
      return F(NAN);

    // It is a good idea to select a proper initial approximation. For  that,
    // extract the binary exponent from the "x" rep; "frexp" is a "constexpr"
    // in C++ >= 23:
    F

  }
  // End namespace "ConstExprMaths"
}
// End namespace SpaceBallistics
