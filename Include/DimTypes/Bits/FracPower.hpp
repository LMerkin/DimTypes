//===========================================================================//
//                       "DimTypes/Bits/FracPower.hpp":                      //
//               Compile-Time Reduction of Fractional Powers                 //
//===========================================================================//
// NB: This header is typically to be included from "Dimensions.h" but can also
// be used stand-alone:
#pragma  once
#include <cstdlib>
#include <cmath>
#include <complex>
#include <cassert>

namespace DimTypes
{
namespace Bits
{
  //=========================================================================//
  // Overloaded / templated functions on representation types:               //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "SqRt":                                                                 //
  //-------------------------------------------------------------------------//
  inline float SqRt (float x)
    { return ::sqrtf(x); }

  inline double SqRt(double x)
    { return ::sqrt (x); }

  inline long double SqRt(long double x)
    { return ::sqrtl(x); }

  template<typename T>
  inline std::complex<T> SqRt(std::complex<T> z)
    { return std::sqrt(z); }

  //-------------------------------------------------------------------------//
  // "CbRt":                                                                 //
  //-------------------------------------------------------------------------//
  inline float CbRt(float x)
    { return ::cbrtf(x); }

  inline double CbRt(double x)
    { return ::cbrt(x);  }

  inline long double CbRt(long double x)
    { return ::cbrtl(x); }

  // NB: There is no built-in "CbRt" for "complex" types, so use the "pow":
  template<typename T>
  inline std::complex<T> CbRt(std::complex<T> x)
    { return std::pow(x, T(1.0)/T(3.0)); }

  //-------------------------------------------------------------------------//
  // "Abs":                                                                  //
  //-------------------------------------------------------------------------//
  inline float Abs(float x)
    { return ::fabsf(x); }

  inline double Abs(double x)
    { return ::fabs(x);  }

  inline long double Abs(long double x)
    { return ::fabsl(x); }

  // NB: For "complex" types, the result of "Abs" is lifted back to "complex":
  template<typename T>
  inline std::complex<T>      Abs(std::complex<T> z)
    { return  std::complex<T>(abs(z)); }

  //-------------------------------------------------------------------------//
  // "Floor":                                                                //
  //-------------------------------------------------------------------------//
  inline float Floor(float x)
    { return ::floorf(x); }

  inline double Floor(double x)
    { return ::floor(x);  }

  inline long double Floor(long double x)
    { return ::floorl(x); }

  // XXX: For "complex" types, "Floor" is applied to both Re and Im parts:
  template<typename T>
  inline std::complex<T>      Floor(std::complex<T>  z)
    { return  std::complex<T>(Floor(z.real()), Floor(z.imag())); }

  //-------------------------------------------------------------------------//
  // "Ceil":                                                                 //
  //-------------------------------------------------------------------------//
  inline float Ceil(float x)
    { return ::ceilf(x); }

  inline double Ceil(double x)
    { return ::ceil(x);  }

  inline long double Ceil(long double x)
    { return ::ceill(x); }

  // XXX: For "complex" types, "Ceil" is applied to both Re and Im parts:
  template<typename T>
  inline std::complex<T>      Ceil(std::complex<T> z)
    { return  std::complex<T>(Ceil(z.real()), Ceil(z.imag())); }

  //-------------------------------------------------------------------------//
  // "Round":                                                                //
  //-------------------------------------------------------------------------//
  inline float Round(float x)
    { return ::roundf(x); }

  inline double Round(double x)
    { return ::round(x);  }

  inline long double Round(long double x)
    { return ::roundl(x); }

  // NB: For "complex" types, "Round" is applied to both Re and Im parts:
  template<typename T>
  inline std::complex<T>      Round(std::complex<T> z)
    { return  std::complex<T>(Round(z.real()), Round(z.imag())); }

  //=========================================================================//
  // Compile-Time GCD and Normalisation Functions:                           //
  //=========================================================================//
  constexpr int GCDrec(unsigned p, unsigned q)
  {
    assert (p <= q);
    return (p == 0) ? q : GCDrec(q % p, p);
  }

  constexpr int GCD(int m, int n)
  {
		// Normalise both "m" and "n" first; GCD will always be > 0:
    unsigned p = unsigned((m >= 0) ? m : (-m));
    unsigned q = unsigned((n >= 0) ? n : (-n));
    return  (p <= q)
             ? GCDrec(p, q)
             : GCDrec(q, p);
  }

  constexpr int NormaliseNumer(int m, int n)
  {
    // NB: GCD is always >= 0. Compile-time error occurs if GCD == 0 (which can
    // happen only if m == n == 0):
    return (n > 0)
            ?   (m / GCD(m, n))
            : - (m / GCD(m, n));
  }

  constexpr int NormaliseDenom(int m, int n)
  {
    int    an = (n >= 0) ? n : (-n);
    return an / GCD(m, n);
  }

  //=========================================================================//
  // Power Templated Functions:                                              //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // Power Expr with a Natural Degree:                                       //
  //-------------------------------------------------------------------------//
  // The Even Positive case:
  template<typename T, int M, bool IsEven = (M % 2 == 0)>
  struct NatPower
  {
    static_assert((M >= 2) && IsEven, "Invalid degree in NatPower");
    constexpr static T res(T base)
    {
      T halfPow = NatPower<T, M/2>::res(base);
      return halfPow * halfPow;
    }
  };

  // The Odd case:
  template<typename T, int M>
  struct NatPower<T, M, false>
  {
    static_assert(M >= 1, "Invalid degree in NatPower");
    constexpr static T res(T base)
    {
      T halfPow = NatPower<T, M/2>::res(base);
      return base * halfPow * halfPow;
    }
  };

  // M==1: This case is not really necessary, but is provided for efficiency:
  template<typename T>
  struct NatPower<T, 1, false>
  {
    constexpr static T res(T base)
      { return base; }
  };

  // M==0: Trivial case:
  template<typename T>
  struct NatPower<T, 0, true>
  {
    constexpr static T res(T)
      { return T(1.0); }  // XXX: We do not check if base==0
  };

  //-------------------------------------------------------------------------//
  // Power Expr with a General Integral Degree:                              //
  //-------------------------------------------------------------------------//
  // Natural case:
  template<typename T, int M, bool IsNatural = (M >= 0)>
  struct IntPower
  {
    static_assert(IsNatural, "IntPower: Case match error");
    constexpr static T res(T base)
      { return NatPower<T, M>::res(base); }
  };

  // The Negative Integral case:
  template<typename T, int M>
  struct IntPower<T, M, false>
  {
    constexpr static T res(T base)
      { return T(1.0) / NatPower<T, -M>::res(base); }
  };

  //-------------------------------------------------------------------------//
  // "SqRtPower":                                                            //
  //-------------------------------------------------------------------------//
  // The degree is a (normalised) Rational; testing the denom "N" for being a
  // multiple of 2 -- if so, rewriting the power expr via "SqRt".
  //
  // Generic case: N != 1 and N is NOT a multiple of 2: use "pow":
  template<typename T, int M, unsigned N, bool IsSqrt = (N % 2 == 0)>
  struct SqRtPower
  {
    static_assert(N >= 3 && !IsSqrt, "SqRtPower: Degree not normalised");
    static T res(T base)
      { return std::pow(base, T(M) / T(N)); }
  };

  // The Integral degree case (N==1):
  template<typename T, int M>
  struct SqRtPower<T, M, 1, false>
  {
    constexpr static T res(T base)
      { return IntPower<T, M>::res(base); }
  };

  // "Sqrt" case: NB: "N" is indeed a multiple of 2:
  template<typename T, int M, unsigned N>
  struct SqRtPower<T, M, N, true>
  {
    static_assert(N >= 2 && N % 2 == 0, "SqRtPower: Case match error");
    constexpr static T res(T base)
      { return SqRtPower<T, M, N/2>::res(SqRt(base)); }
  };

  //-------------------------------------------------------------------------//
  // "CbRtPower":                                                            //
  //-------------------------------------------------------------------------//
  // Similar to "SqRtPower"; testing "N" for being a multiple of 3,  and if so,
  // rewriting the power expr via "CbRt". However, "CbRt" is only available for
  // real types.
  //
  // Generic case: N is NOT a multiple of 3, or inappropritate type "T":   Fall
  // back to "SqRtPower":
  template<typename T, int M, unsigned N, bool IsCbrt = (N % 3 == 0)>
  struct CbRtPower
  {
    static T res(T base)
      { return SqRtPower<T, M, N>::res(base); }
  };

  // The following types allow us to use "CbRt":
  template<int M, unsigned N>
  struct CbRtPower<float,  M, N, true>
  {
    static_assert(N >= 3 && N % 3 == 0, "CbRtPower: Case match error");
    static float res(float base)
      { return CbRtPower<float, M, N/3>::res(CbRt(base)); }
  };

  template<int M, unsigned N>
  struct CbRtPower<double, M, N, true>
  {
    static_assert(N >= 3 && N % 3 == 0, "CbRtPower: Case match error");
    static double res(double base)
      { return CbRtPower<double, M, N/3>::res(CbRt(base)); }
  };

  template<int M, unsigned N>
  struct CbRtPower<long double, M, N, true>
  {
    static_assert(N >= 3 && N % 3 == 0, "CbRtPower: Case match error");
    static long double res(long double base)
      { return CbRtPower<long double, M, N/3>::res(CbRt(base)); }
  };

  //-------------------------------------------------------------------------//
  // Power Expr with a General Fractional Degree:                            //
  //-------------------------------------------------------------------------//
  // NB: This is a template function, not a struct -- it has no partial specs.
  // First attempt "CbRtPower", then it gets reduced further:
  template<int M, unsigned N, typename T>
  inline T FracPower(T base)
  {
    static_assert(N != 0, "Zero denom in fractional degree");
    return CbRtPower<T, NormaliseNumer(M,N), NormaliseDenom(M,N)>::res(base);
  }
}
// End namespace Bits
}
// End namespace DimTypes
