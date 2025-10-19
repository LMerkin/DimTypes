// vim:ts=2:et
//===========================================================================//
//                        "DimTypes/Bits/Encodings.hpp":                     //
//         Encoding/Decoding of Dimension Exponents and Unit Vectors         //
//===========================================================================//
#pragma once
#include "CEMaths.hpp"
#include "Macros.h"
#include <complex>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstdint>

namespace DimTypes
{
namespace Bits
{
  //=========================================================================//
  // "Encodings": Zp Representation of Rationals, and more:                  //
  //=========================================================================//
  // Dimension Exponents are monomials over Fundamental Dims with Rational pows.
  // They will be represented as "uint64_t" in which bit flds of "PBits" will be
  // allocated for each dimension. Bit fields encode a rational powers in the Zp
  // format where p = "PMod", which is the largest prime fitting in the "PBits"
  // format, ie (PMod < 2**PBits):
  // We allow    "MaxDims" to be:   7,   8,   9;
  // the corresp "PBits"      is:   9,   8,   7;
  // the corresp "PMod"       is: 509, 251, 127;
  // All Unit enums must be encodeable in the "PBits" format; but this is not a
  // real constraint at the moment (can create up to "PMod" different units per
  // dim), so the number of units is not checked:
  //
  constexpr inline unsigned DefMaxDims = 8;

  template<typename RepT, unsigned MaxDims>
  class Encodings
  {
  public:
    static_assert(MaxDims  == 7 || MaxDims == 8  || MaxDims == 9);

    constexpr static unsigned PBits =
      (MaxDims == 7) ? 9   : (MaxDims == 8) ? 8   : 7;

    constexpr static uint64_t PMod  =
      (MaxDims == 7) ? 509 : (MaxDims == 8) ? 251 : 127;

    constexpr static int      IPMod = int(PMod);

    constexpr static uint64_t PMask = (1UL << PBits) - 1;
    static_assert(PMod <= PMask, "Insufficient bits for PMod");

    static_assert(PBits * MaxDims <= sizeof(uint64_t) * 8,
                  "Exp Vector does not fit into uint64_t");

    // Exponents for Fundamental Dimensions: Each exponent is 1 for the corresp
    // Dim:
    constexpr static uint64_t DimExp(unsigned dim)
    {
      assert(dim <  MaxDims);
      return 1UL << (dim * PBits);
    }

    //=======================================================================//
    // Compile-Time GCD and Normalisation Functions:                         //
    //=======================================================================//
  private:
    constexpr static unsigned GCDrec(unsigned p, unsigned q)
    {
      assert (p <= q);
      return (p == 0) ? q : GCDrec(q % p, p);
    }

  public:
    constexpr static unsigned GCD(int m, int n)
    {
      // Normalise both "m" and "n" first; GCD will always be > 0:
      unsigned p = unsigned((m >= 0) ? m : (-m));
      unsigned q = unsigned((n >= 0) ? n : (-n));
      return  (p <= q)
               ? GCDrec(p, q)
               : GCDrec(q, p);
    }

  private:
    constexpr static std::pair<int, unsigned> NormaliseFrac(int m, int n)
    {
      assert(n != 0);
      // NB: GCD is always >= 0. Compile-time error occurs if GCD == 0 (which
      // can happen only if m == n == 0):
      int gcd = int(GCD(m, n));
      int m1  = (n > 0 ? m : (-m)) / gcd;
      int n1  = (n > 0 ? n : (-n)) / gcd;
      assert(n1 > 0);
      return std::make_pair(m1, unsigned(n1));
    }

  public:
    // Normalising a quantity (positive or negative) modulo PMod.
    // The result is always in [0 .. PMod-1]:
    constexpr static unsigned Normalise(int x)
    {
      int    res = (x >= 0) ? (x % IPMod) : (x % IPMod + IPMod);
      assert(res >= 0);
      return unsigned(res);
    }

  private:
    // "InverseModP":
    // Zp inverse using the Extended GCD algorithm. Returns the coeff "c" s.t.
    // GCD(x,P) = 1 = c*x + d*P, c >= 0, where P = IPMod
    // Pre-condition: 0 <= x && x < y ;
    // Result:    in [0 .. PMod-1]    :
    //
    constexpr static unsigned InverseModP(int n)
    {
      if (n % IPMod == 0)
        throw "DimTypes::Bits::InverseModP: ERROR: UnInvertible arg";

      int x = int(Normalise(n));
      int a = 1;
      int b = 0;
      int y = IPMod;
      int c = 0;
      int d = 1;
      // Check failure will result in type error at compile time:
      if (!(0 <= x && x < y))
        throw "DimTypes::Bits::InverseModP: LOGIC ERROR (1)";

      while (x != 0)
      {
        int q = y / x;
        int r = y % x;
        y = x;
        x = r;
        if (!(0 <= x && x < y))  // As above
          throw "DimTypes::Bits::InverseModP: LOGIC ERROR (2)";
        int a1 = c - q * a;
        int b1 = d - q * b;
        c = a;
        d = b;
        a = a1;
        b = b1;
      }
      return Normalise(c);
    }

  public:
    //=======================================================================//
    // Compile-Time Monomial Operations on Dimension Exponents:              //
    //=======================================================================//
    // Addition, Subtraction, Multiplication by Rationals:
    //
    constexpr static uint64_t GetFld(uint64_t From, unsigned dim)
    {
      // Move the selected bit field to the right and zero-out all other bits:
      assert (dim  < MaxDims);
      return (From >> (dim * PBits)) & PMask;
    }

  private:
    constexpr static uint64_t PutFld(uint64_t From, unsigned dim)
    {
      // Zero out the upper bits and move lower ones to the left into the
      // required position:
      assert (dim  < MaxDims);
      return (From & PMask) << (dim * PBits);
    }

  public:
    constexpr static uint64_t AddExp(uint64_t E, uint64_t F)
    {
      // Adding up bit flds of the exponents "E" and "F" modulo "PMod":
      if (E == 0)
        return F;
      if (F == 0)
        return E;

      uint64_t res = 0UL;
      for (unsigned dim = 0;  dim < MaxDims;  ++dim)
        res |= PutFld((GetFld(E, dim) + GetFld(F, dim)) % PMod, dim);
      return res;
    }

    constexpr static uint64_t SubExp(uint64_t E, uint64_t F)
    {
      // Subtracting bit flds of the exponents "E" and "F" modulo "PMod":
      if (F == 0)
        return E;

      uint64_t res = 0UL;
      for (unsigned dim = 0;  dim < MaxDims;  ++dim)
        // NB: First, add "PMod" to the LHS to avoid negative vals:
        res |= PutFld(((PMod + GetFld(E, dim)) - GetFld(F, dim)) % PMod, dim);
      return res;
    }

    constexpr static uint64_t MultExp(uint64_t E, int m)
    {
      // Multiplying bit flds of the exponent "E" by "m" modulo "PMod":
      if (m == 1)
        return E;

      uint64_t res = 0UL;
      for (unsigned dim = 0;  dim < MaxDims;  ++dim)
        res |= PutFld((GetFld(E, dim) * uint64_t(Normalise(m))) % PMod, dim);
      return res;
    }

    constexpr static uint64_t DivExp(uint64_t E, unsigned n)
    {
      // Dividing bit flds of the exponent "E" by "n" modulo "PMod":
      if (n == 0)
        throw "DimTypes::Bits::DivExp: LOGIC ERROR";
      if (n == 1)
        return E;

      uint64_t res = 0UL;
      for (unsigned dim = 0;  dim < MaxDims;  ++dim)
        res |=
          PutFld((GetFld(E, dim) * uint64_t(InverseModP(int(n)))) % PMod, dim);
      return res;
    }

    //=======================================================================//
    // "MaxHeight", "FindMaxHeight":                                         //
    //=======================================================================//
  private:
    // It turns out that the maximum height of a fraction respresentable in the
    // mod IPMod format (without collisions) is much less than IPMod:
    //
    // Prevent bogus warnings in CLang (this function is invoked at compile-
    // time, so the compiler would detect any erroneous behavior anyway):
    constexpr static unsigned FindMaxHeight()
    {
      PREVENT_BOGUS_CLANG_WARNINGS_B
      bool taken[IPMod];
      for (int i = 0; i < IPMod; ++i)
        taken[i] = false;

      for (int height = 2;  height < IPMod;  ++height)
      for (int denom  = 1;  denom  < height; ++denom)
      {
        unsigned invDenom = InverseModP(denom);  // XXX: Repeated computation
        int      numerP   = height - denom;      // Pos numer

        // However, consider only irreducible fractions:
        if (GCD(numerP, denom) != 1)
          continue;

        int      numerC   = IPMod  - numerP;     // Complement of Neg numer
        if (!(numerP > 0 && numerC > 0))
          throw "DimTypes::Bits::FindMaxHeight: LOGIC ERROR";

        unsigned repP     = (unsigned(numerP) * invDenom) % PMod;
        unsigned repC     = (unsigned(numerC) * invDenom) % PMod;

        if (taken[repP] || taken[repC])
          return unsigned(height-1); // Clash encountered at "height"!

        taken[repP] = true;
        taken[repC] = true;
      }
      return PMod-1;       // We will not really get here...
      PREVENT_BOGUS_CLANG_WARNINGS_E
    }

  public:
    // So:
    constexpr static unsigned MaxHeight = FindMaxHeight();

    //=======================================================================//
    // Extraction of Numer and Denom from Zp-Encodings:                      //
    //=======================================================================//
    // Find the Numer and Denom of the minimal total height  corresponding  to
    // the given Zp representation. XXX: Currently this is done just by direct
    // search:
    //
    constexpr static std::pair<int, unsigned> GetNumerAndDenom(uint64_t rep)
    {
      // Numer is +-(height - denom):
      if (rep == 0UL)
        return std::make_pair(0, 1);

      // Generic case: Traverse all heights. Typically, the result is nearby:
      // height == |numer| + denom, where denom >= 1 and numer != 0:
      for (int height = 2; height < int(MaxHeight); ++height)
      {
        for (int denom = 1; denom < height; ++denom)
        {
          unsigned invDenom = InverseModP(denom);
          int      numerP   = height - denom;  // Positive numer
          int      numerC   = IPMod  - numerP; // Complement of Negative numer

          if (!(numerP > 0 && numerC > 0))
            throw "DimTypes::Bits::GetNumerAndDenom: LOGIC ERROR";

          if ((unsigned(numerP) * invDenom) % PMod == rep)
            return std::make_pair(numerP,   unsigned(denom));
          if ((unsigned(numerC) * invDenom) % PMod == rep)
            return std::make_pair(-numerP,  unsigned(denom));
        }
      }
      // XXX: Can we get here without a success. Throwing an exception will
      // result in a type error at compile time:
      throw "ERROR: DimTypes::Bits::GetNumerAndDenom: Rep Not Matched";
    }

    //=======================================================================//
    // Compile-Time Operations on Unit Vectors:                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Setting / Creation of Unit Encodings:                                 //
    //-----------------------------------------------------------------------//
    constexpr static uint64_t SetUnit(uint64_t U, unsigned dim, unsigned unit)
    {
      assert(dim < MaxDims);
      // Clear the old bits at "dim" and set those from "unit" there.
      // XXX: we must always have unit <= PMask,  but instead of checking this
      // pre-cond, we simply apply "PMask" to "unit". In all normal cases this
      // should be fine:
      return (U   & ~(PMask  << (dim * PBits))) |
             ((unit & PMask) << (dim * PBits));
    }

    constexpr static uint64_t MkUnit(unsigned dim, unsigned unit)
    {
      assert(dim < MaxDims);
      return SetUnit(0UL, dim, unit);
    }

    //-----------------------------------------------------------------------//
    // Unification of Units:                                                 //
    //-----------------------------------------------------------------------//
    // See the implementation for the exact semantics. Unifies the units in the
    // (E,U) and (F,V) operands of "*" or "/":
    //
    constexpr static uint64_t UnifyUnits
      (uint64_t E, uint64_t F, uint64_t U, uint64_t V)
    {
      uint64_t res = 0UL;
      for (unsigned dim = 0; dim < MaxDims; ++dim)
      {
        uint64_t e = GetFld(E, dim);
        uint64_t f = GetFld(F, dim);
        uint64_t u = GetFld(U, dim);
        uint64_t v = GetFld(V, dim);

        // Unified Units:
        uint64_t unified =
          (e == 0UL)
          ? // "e" is dimension-less:
            ((f == 0UL)
             ? // "f" is dimension-less as well, so both units are reset:
               0L
             : // Use the unit of "f":
               v
            )
          : // "e" is non-trivial:
            ((f == 0UL)
             ? // "f" is dimension-less, so use the unit of "e":
               u
             : // Both "e" and "f" are non-trivial, so their units must be
               // same:
               (u == v)
               ? u
               : throw "ERROR: DimTypes::Bits::UnifyUnits: Unification Failed"
            );
        // Put the unified units into the "res":
        res |= PutFld(unified, dim);
      }
      return res;
    }

    // "UnitsOK": Similar for "UnifyUnits" but for one Exponent Vector.
    // Returns a boolean used in external static assertions:
    //
    constexpr static bool UnitsOK(uint64_t E, uint64_t U, uint64_t V)
    {
      for (unsigned dim = 0; dim <  MaxDims; ++dim)
        if (GetFld(E,dim) != 0UL && GetFld(U,dim) != GetFld(V,dim))
          // Failed units unification: Exp != 0, and units differ:
          return false;
      // If we got here:
      return true;
    }

    //-----------------------------------------------------------------------//
    // Clean-Up (Re-Setting to the default 0) of unused units:               //
    //-----------------------------------------------------------------------//
    constexpr static uint64_t CleanUpUnits(uint64_t E, uint64_t U)
    {
      uint64_t res = 0UL;
      for (unsigned dim = 0; dim < MaxDims; ++dim)
      {
        uint64_t u =
          (GetFld(E, dim) != 0)
          ? // This dim's Exp is non-trivial, so the Units are indeed required:
            GetFld(U, dim)
          : 0UL;
        res |= PutFld(u, dim);
      }
      return res;
    }

    //=======================================================================//
    // Power Templated Functions:                                            //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Power Expr with a Natural Degree:                                     //
    //-----------------------------------------------------------------------//
    template<int M>
    constexpr static RepT IntPow(RepT a_x)
    {
      if constexpr(M < 0)
        return RepT(1.0) / IntPow<-M>(a_x);

      if constexpr(M == 0)
        return RepT(1.0);

      if constexpr(M == 1)
        return a_x;

      RepT halfPow  = IntPow<M/2>(a_x);
      RepT halfPow2 = halfPow * halfPow;
      if constexpr(M % 2 == 1)
        return halfPow2 * a_x;
      else
        return halfPow2;
    }

  private:
    //-----------------------------------------------------------------------//
    // "Only2and3":                                                          //
    //-----------------------------------------------------------------------//
    constexpr static bool Only2and3(unsigned a_n)
    {
      return
        (a_n == 0 || a_n == 1)
        ? true :
        (a_n % 2 == 0)
        ? Only2and3(a_n / 2) :
        (a_n % 3 == 0)
        ? Only2and3(a_n / 3)
        : false;
    }

    //-----------------------------------------------------------------------//
    // "FracPow23":                                                          //
    //-----------------------------------------------------------------------//
    // "N" is assumed to consist of 2 and 3 multiples only, so we can use the
    // square and cubic roots:
    //
    template<int M, unsigned N>
    constexpr static RepT FracPow23(RepT a_x)
    {
      static_assert(M != 0 && N != 0 && Only2and3(N),
                    "FracPow23: Invalid Frac");

      if constexpr(N == 1)
        return IntPow<M>(a_x);
      else
      if constexpr(N % 2 == 0)
        return FracPow23<M, N/2>(CEMaths::SqRt<RepT>(a_x));
      else
      {
        static_assert(N % 3 == 0, "FracPow23: N != Mults(2,3)");
        return FracPow23<M, N/3>(CEMaths::CbRt<RepT>(a_x));
      }
    }

  public:
    //-----------------------------------------------------------------------//
    // Power Expr with a General Fractional Degree:                          //
    //-----------------------------------------------------------------------//
    template<int M, unsigned N>
    constexpr static RepT FracPow(RepT a_x)
    {
      constexpr auto     nFrac = NormaliseFrac(M, N);
      constexpr int      M1    = nFrac.first;
      constexpr unsigned N1    = nFrac.second;

      // Integral power?
      if constexpr(M1 == 0)
        return RepT(1.0);
      if constexpr(N1 == 1)
        return IntPow<M1>(a_x);

      // If the power is indeed fractional, we check whether "N1" only consists
      // of 2 and 3 multiples,  in which case, simplifications  based on "SqRt"
      // and "CbRt" may be available:
      //
      if constexpr(Only2and3(N1))
        return FracPow23<M1, N1>(a_x);
      else
      if constexpr(CEMaths::IsComplex<RepT>)
      {
        // If "RepT" is "complex", the type of the degree is still the underly-
        // ing real one:
        using P = typename  RepT::value_type;
        return CEMaths::Pow<RepT>(a_x, P(M1)/P(N1));
      }
      else
        // Otherwise, use the generic real "Pow":
        return CEMaths::Pow<RepT>(a_x, RepT(M1) / RepT(N1));
    }

    //=======================================================================//
    // Misc:                                                                 //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // "PutMagnitude":                                                       //
    //-----------------------------------------------------------------------//
    // Outputs a Real or Complex magnitude to the given buffer:
    //
    // Generic Case:  A Real Value: convert it to "double":
    static char* PutMagnitude(char* a_buff, int a_n, RepT const& a_val)
    {
      PREVENT_BOGUS_CLANG_WARNINGS_B
      assert(a_buff != nullptr && a_n > 0);

      if constexpr(CEMaths::IsComplex<RepT>)
      {
        // Special Case:  A Complex Value:
        double  magRe  = double(a_val.real());
        double  magIm  = double(a_val.imag());

        return a_buff +  snprintf(a_buff, size_t(a_n), "(%.16e %c %.16e * I)",
                         magRe, (magIm < 0.0) ? '-'      : '+',
                                (magIm < 0.0) ? (-magIm) : magIm);
      }
      else
        // Generic Case:  A Real Value:
        return a_buff + snprintf(a_buff, size_t(a_n), "%.16e", double(a_val));

      PREVENT_BOGUS_CLANG_WARNINGS_E
    }
  };
  // End "Encodings" class
}
// End namespace Bits
}
// End namespace DimTypes
