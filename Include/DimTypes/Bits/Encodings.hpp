// vim:ts=2:et
//===========================================================================//
//                        "DimTypes/Bits/Encodings.hpp":                     //
//         Encoding/Decoding of Dimension Exponents and Unit Vectors         //
//===========================================================================//
#pragma once
#include <complex>

namespace DimTypes
{
namespace Bits
{
  //=========================================================================//
  // Zp Representation of Rationals:                                         //
  //=========================================================================//
  // Dimension Exponents are monomials over Fundamental Dims with Rational pows.
  // They will be represented as "unsigned long" integers in which bit flds  of
  // "PBits" will be allocated for each dimension. Bit fields encode a rational
  // powers in the Zp format where p = PMod below. Thus, 7 dims must fit into a
  // 64-bit format -- so PBits = 9 and PMod = 509 (the largest prime which fits
  // into 9 bits, ie is smaller than 2^9=512):
  // All Unit enums must be encodeable in the "PBits" format; but this is not a
  // real constraint at the moment (can create up to "PMod" different units per
  // dim), so the number of units is not checked.
  //
  using                   ULong = unsigned long;  // To avoid warnings in CLang
  constexpr unsigned long PMod  = 509;
  constexpr int IPMod           = int(PMod);

  constexpr unsigned      NDims = 7;
  constexpr unsigned      PBits = 9;
  constexpr unsigned long PMask = (1UL << PBits) - 1;
  static_assert(PMod <= PMask, "Insufficient bits for PMod");

  static_assert(PBits * NDims <= sizeof(unsigned long) * 8,
                "Exp Vector does not fit into unsigned long");

  // Exponents for Fundamental Dimensions: Each exponent is 1 for the corresp
  // Dim:
  constexpr unsigned long DimExp(unsigned dim)
    { return 1UL << (dim * PBits); }

  //=========================================================================//
  // Compile-Time Monomial Operations on Dimension Exponents:                //
  //=========================================================================//
  // Addition, Subtraction, Multiplication by Rationals:
  // XXX: in the functions below, it would be better to check that "dim" is in
  // the valid range (0..NDims-1), but generating proper compile-time err msgs
  // is difficult, so we don't do it at the moment:
  //
  constexpr unsigned long GetFld(unsigned long From, unsigned dim)
  {
    // Move the selected bit field to the right and zero-out all other bits:
    return (From >> (dim * PBits)) & PMask;
  }

  constexpr unsigned long PutFld(unsigned long From, unsigned dim)
  {
    // Zero out the upper bits and move lower ones to the left into the
    // required position:
    return (From & PMask) << (dim * PBits);
  }

  // Normalising a quantity (positive or negative) modulo PMod.
  // The result is always in [0 .. PMod-1]:
  constexpr unsigned Normalise(int x)
    { return (x >= 0) ? (x % IPMod) : (x % IPMod + IPMod); }

  // Zp inverse using the Extended GCD algorithm.
  // The following function returns 
  // coeffs "c" such that
  // GCD(x,y) = 1 = c*x + d*y, c >= 0, where y = IPMod
  // Pre-condition: 0 <= x && x < y  :
  // Result:     in [0 .. PMod-1]:
  //
  constexpr unsigned InverseModP(int n)
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

  constexpr unsigned long AddExp(unsigned long E, unsigned long F)
  {
    // Adding up bit flds of the exponents "E" and "F" modulo "PMod":
    if (E == 0)
      return F;
    if (F == 0)
      return E;

    unsigned long res = 0UL;
    for (unsigned dim = 0;  dim < NDims;  ++dim)
      res |= PutFld((GetFld(E, dim) + GetFld(F, dim)) % PMod, dim);
    return res;
  }

  constexpr unsigned long SubExp(unsigned long E, unsigned long F)
  {
    // Subtracting bit flds of the exponents "E" and "F" modulo "PMod":
    if (F == 0)
      return E;

    unsigned long res = 0UL;
    for (unsigned dim = 0;  dim < NDims;  ++dim)
      // NB: First, add "PMod" to the LHS to avoid negative vals:
      res |= PutFld(((PMod + GetFld(E, dim)) - GetFld(F, dim)) % PMod, dim);
    return res;
  }

  constexpr unsigned long MultExp(unsigned long E, int m)
  {
    // Multiplying bit flds of the exponent "E" by "m" modulo "PMod":
    if (m == 1)
      return E;

    unsigned long res = 0UL;
    for (unsigned dim = 0;  dim < NDims;  ++dim)
      res |= PutFld((GetFld(E, dim) * ULong(Normalise(m))) % PMod, dim);
    return res;
  }

  constexpr unsigned long DivExp(unsigned long E, unsigned n)
  {
    // Dividing bit flds of the exponent "E" by "n" modulo "PMod":
    if (n == 0)
      throw "DimTypes::Bits::DivExp: LOGIC ERROR";
    if (n == 1)
      return E;

    unsigned long res = 0UL;
    for (unsigned dim = 0;  dim < NDims;  ++dim)
      res |= PutFld((GetFld(E, dim) * ULong(InverseModP(n))) % PMod, dim);
    return res;
  }

  //=========================================================================//
  // "MaxHeight", "FindMaxHeight":                                           //
  //=========================================================================//
  // It turns out that the maximum height of a fraction respresentable in the
  // mod IPMod format (without collisions) is much less than IPMod:
  //
  constexpr unsigned FindMaxHeight()
  {
    bool taken[IPMod];
    for (int i = 0; i < IPMod; ++i)
      taken[i] = false;

    for (int height = 2;  height < IPMod;  ++height)
    for (int denom  = 1;  denom  < height; ++denom)
    {
      unsigned invDenom = InverseModP(denom);      // XXX: Repeated computations
      int      numerP   = height - denom;          // Pos numer

      // However, consider only irreducible fractions:
      if (GCD(numerP, denom) != 1)
        continue;

      int      numerC   = IPMod  - numerP;         // Complement of Neg numer
        if (!(numerP > 0 && numerC > 0))
          throw "DimTypes::Bits::FindMaxHeight: LOGIC ERROR";

      unsigned repP     = (unsigned(numerP) * invDenom) % PMod;
      unsigned repC     = (unsigned(numerC) * invDenom) % PMod;

      if (taken[repP] || taken[repC])
        return height-1; // Clash encountered at "height"!

      taken[repP] = true;
      taken[repC] = true;
    }
    return PMod-1;       // We will not really get here...
  }

  // So:
  constexpr unsigned MaxHeight = FindMaxHeight();

  //=========================================================================//
  // Extraction of Numer and Denom from Zp-Encodings:                        //
  //=========================================================================//
  // Find the Numer and Denom of the minimal total height corresponding to the
  // given Zp representation. Currently this is done just by direct search:
  //
  constexpr std::pair<int, unsigned> GetNumerAndDenom(unsigned long rep)
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
        int invDenom = InverseModP(unsigned(denom));
        int numerP   = height - denom;          // Positive numer
        int numerC   = IPMod  - numerP;         // Complement of Negative numer

        if (!(numerP > 0 && numerC > 0))
          throw "DimTypes::Bits::GetNumerAndDenom: LOGIC ERROR";

        if (unsigned(numerP * invDenom) % IPMod == rep)
          return std::make_pair(numerP,   unsigned(denom));
        if (unsigned(numerC * invDenom) % IPMod == rep)
          return std::make_pair(-numerP,  unsigned(denom));
      }
    }
    // XXX: Can we get here without a success. Throwing an exception will
    // result in a type error at compile time:
    throw "ERROR: DimTypes::Bits::GetNumerAndDenom: Rep Not Matched";
  }

  //=========================================================================//
  // Compile-Time Operations on Unit Vectors:                                //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // Setting / Creation of Unit Encodings:                                   //
  //-------------------------------------------------------------------------//
  constexpr unsigned long SetUnit(unsigned long U, unsigned dim, unsigned unit)
  {
    // Clear the old bits at "dim" and set those from "unit" there.
    // XXX: we must always have unit <= PMask,  but instead of checking this
    // pre-cond, we simply apply "PMask" to "unit". In all normal cases this
    // should be fine:
    return (U   & ~(PMask  << (dim * PBits))) |
           ((unit & PMask) << (dim * PBits));
  }

  constexpr unsigned long MkUnit(unsigned dim, unsigned unit)
    { return SetUnit(0UL, dim, unit); }

  //-------------------------------------------------------------------------//
  // Unification of Units:                                                   //
  //-------------------------------------------------------------------------//
  // See the implementation for the exact semantics. Unifies the units in (E,U)
  // and (F,V) operands of "*" or "/":
  //
  constexpr  unsigned long UnifyUnits
    (unsigned long E, unsigned long F, unsigned long U, unsigned long V)
  {
    unsigned long res = 0UL;
    for (unsigned dim = 0; dim < NDims; ++dim)
    {
      unsigned long e = GetFld(E, dim);
      unsigned long f = GetFld(F, dim);
      unsigned long u = GetFld(U, dim);
      unsigned long v = GetFld(V, dim);

      // Unified Units:
      unsigned long unified =
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
           : // Both "e" and "f" are non-trivial, so their units must be same:
             (u == v)
             ? u
             : throw "ERROR: DimTypes::Bits::UnifyUnits: Unification Failed"
          );
      // Put the unified units into the "res":
      res |= PutFld(unified, dim);
    }
    return res;
  }

  // "UnitsOK": Similar for "UnifyUnits" but for one Exponent Vector. Returns a
  // boolean used in external static assertions:
  //
  constexpr bool UnitsOK(unsigned long E, unsigned long U, unsigned long V)
  {
    for (unsigned dim = 0; dim < NDims; ++dim)
      if (GetFld(E,dim) != 0UL && GetFld(U,dim) != GetFld(V,dim))
        // Failed units unification: Exp != 0, and units differ:
        return false;
    // If we got here:
    return true;
  }

  //-------------------------------------------------------------------------//
  // Clean-Up (Re-Setting to the default 0) of unused units:                 //
  //-------------------------------------------------------------------------//
  constexpr unsigned long CleanUpUnits(unsigned long E, unsigned long U)
  {
    unsigned long res = 0UL;
    for (unsigned dim = 0; dim < NDims; ++dim)
    {
      unsigned long u =
        (GetFld(E, dim) != 0)
        ? // This dim's Exp is non-trivial, so the Units are indeed required:
          GetFld(U, dim)
        : 0UL;
      res |= PutFld(u, dim);
    }
    return res;
  }

  //=========================================================================//
  // Misc:                                                                   //
  //=========================================================================//
  // Whether a given type "T" is a complex type:
  //
  template<typename T>
  constexpr bool IsComplex = false;

  template<typename T>
  constexpr bool IsComplex<std::complex<T>> = true;

  /*
  //=========================================================================//
  // "Put":                                                                  //
  //=========================================================================//
  template<unsigned long E>
  char* OutputUnits();

  template<int Numer, bool IsPos = (Numer > 0) >
  struct UnitStrInt
  {
    // Positive integral exponent:
    static int put(char* buff, char const* unit)
      { return sprintf(buff, " %s^%d", unit, Numer); }
  };

  template<int Numer>
  struct UnitStrInt<Numer,false>
  {
    // Negative integral exponent:
    static int put(char* buff, char const* unit)
      { return sprintf(buff, " %s^(%d)", unit, Numer); }
  };

  template<>
  struct UnitStrInt<1,true>
  {
    // Exponent is 1:
    static int put(char* buff, char const* unit)
      { return sprintf(buff, " %s", unit); }
  };

  template<int Numer, unsigned Denom>
  struct UnitStr
  {
    // General rational case:
    static int put(char* buff, char const*  unit)
      { return sprintf(buff, " %s^(%d/%d)", unit, Numer, Denom); }
  };

  template<int Numer>
  struct UnitStr<Numer,1>
  {
    static int put(char* buff, char const*  unit)
      { return UnitStrInt<Numer>::put(buff, unit); }
  };

  template<>
  struct UnitStr<0,1>
  {
    constexpr static int put(char*, char const*)
      { return 0; }  // Nothing to output
  };
  */
}
// End namespace Bits
}
// End namespace DimTypes
