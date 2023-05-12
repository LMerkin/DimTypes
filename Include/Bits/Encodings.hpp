// vim:ts=2
//===========================================================================//
//                             "Bits/Encodings.hpp":                         //
//         Encoding/Decoding of Dimension Exponents and Unit Vectors         //
//===========================================================================//
#pragma  once
#include <cassert>

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
  // powers in the Zp format where p = PMod below. Thus, 9 dims must fit into a
  // 64-bit format -- so PBits = 7 and PMod = 127 (the largest prime which fits
  // in 7 bits):
  // All Unit enums must be encodeable in the "PBits" format; but this is not a
  // real constraint at the moment (can create up to "PMod" different units per
  // dim), so the number of units is not checked.
  //
  using                   ULong = unsigned long;  // To avoid warnings in CLang
  constexpr unsigned long PMod  = 127;
  constexpr int IPMod           = int(PMod);

  constexpr unsigned      NDims = 9;
  constexpr unsigned      PBits = 7;
  constexpr unsigned long PMask = (1UL << PBits) - 1;
  static_assert(PMod <= PMask, "Insufficient bits for PMod");

  static_assert(PBits * NDims <= sizeof(unsigned long) * 8,
                "Exp Vector does not fit into unsigned long");

  // It turns out that the maximum height of a fraction respresentable in the
  // mod IPMod format is much less than IPMod.   XXX: the following fugure is
  // obtained empirically; it will change when "PMod" changes. For Heights in
  // 1..13, there are no clashes; for 14..15, partial clashes; from 16 on, on-
  // clashes (so heights >= 16 are not used):
  constexpr int MaxHeight = 15;

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
  constexpr unsigned long GetFld    (unsigned long from, unsigned dim)
  {
    // Move the selected bit field to the right and zero-out all other bits:
    return (from >> (dim * PBits)) & PMask;
  }

  constexpr unsigned long PutFld    (unsigned long from, unsigned dim)
  {
    // Zero out the upper bits and move lower ones to the required position:
    return (from & PMask) << (dim * PBits);
  }

  constexpr unsigned long ZeroOutFld(unsigned long from, unsigned dim)
  {
    // All other flds are preserved:
    return from & ~(PMask << (dim * PBits));
  }

  constexpr unsigned long AddFldsModP
    (unsigned long left, unsigned long right, unsigned dim)
  {
    return PutFld((GetFld(left, dim) + GetFld(right, dim)) % PMod, dim);
  }

  constexpr unsigned long SubFldsModP
    (unsigned long left, unsigned long right, unsigned dim)
  {
    // NB: "PMod" is added to make sure no negative values occur; it is then
    // removed by the "%" operator anyway:
    return
    PutFld(((PMod + GetFld(left, dim)) - GetFld(right, dim)) % PMod, dim);
  }

  // Normalising a quantity (positive or negative) modulo PMod.
  // The result is always in [0 .. PMod-1]:
  constexpr int Normalise(int x)
    { return (x >= 0) ? (x % IPMod) : (x % IPMod + IPMod); }

  constexpr unsigned long MultFldsModP(unsigned long left, int m, unsigned dim)
  {
    // NB: "m" may be negative, so first reduce it to a non-negative value
    // modulo "PMod":
    return PutFld((GetFld(left, dim) * ULong(Normalise(m))) % PMod, dim);
  }

  // Zp inverse using the Extended GCD algorithm.
  // The following function returns 
  // coeffs "c" such that
  // GCD(x,y) = 1 = c*x + d*y, c >= 0.
  // Pre-condition: 0 <= x && x < y  :
  constexpr int InverseModP(int x, int a, int b, int y, int c, int d)
  {
    return
      (x == 0)
    ? Normalise(c)
    : InverseModP(y % x, c - (y/x) * a, d - (y/x) * b, x, a, b);
  }

  constexpr int InverseModP(int n)
    { return InverseModP(Normalise(n), 1, 0, IPMod, 0, 1); }

  constexpr unsigned long DivFldsModP(unsigned long left, int n, unsigned dim)
    { return PutFld((GetFld(left, dim) * ULong(InverseModP(n))) % PMod, dim); }

  constexpr unsigned long AddExpRec
    (unsigned long E, unsigned long F, unsigned long curr, unsigned dim)
  {
    return
      (dim >= NDims)
    ? curr
    : AddExpRec(E, F, curr | AddFldsModP(E, F, dim), dim + 1);
  }

  constexpr unsigned long AddExp(unsigned long E, unsigned long F)
  {
    // Adding up bit flds of the exponents "E" and "F" modulo "PMod":
    return
      (E == 0UL)  // First, trivial cases for efficiency:
    ? F
    : (F == 0UL)
    ? E
    : AddExpRec(E, F, 0UL, 0);
  }

  constexpr unsigned long SubExpRec
    (unsigned long E, unsigned long F, unsigned long curr, unsigned dim)
  {
    return
      (dim >= NDims)
    ? curr
    : SubExpRec(E, F, curr | SubFldsModP(E, F, dim), dim + 1);
  }

  constexpr unsigned long SubExp(unsigned long E, unsigned long F)
  {
    // Subtracting bit flds of the exponents "E" and "F" modulo "PMod":
    return
      (F == 0UL)  // First, trivial case for efficiency:
    ? E
    : SubExpRec(E, F, 0UL, 0);
  }

  constexpr unsigned long MultExpRec
    (unsigned long E, int m, unsigned long curr, unsigned dim)
  {
    return
      (dim >= NDims)
    ? curr
    : MultExpRec(E, m, curr | MultFldsModP(E, m, dim), dim + 1);
  }

  constexpr unsigned long MultExp(unsigned long E, int m)
  {
    // Multiplying bit flds of the exponent "E" by "m" modulo "PMod":
    return
      (m == 0)  // First, trivial cases for efficiency:
    ? 0UL
    : (m == 1)
    ? E
    : MultExpRec(E, m, 0UL, 0);
  }

  constexpr unsigned long DivExpRec
    (unsigned long E, int n, unsigned long curr, unsigned dim)
  {
    return
      (dim >= NDims)
    ? curr
    : DivExpRec(E, n, curr | DivFldsModP(E, n, dim), dim + 1);
  }

  constexpr unsigned long DivExp(unsigned long E, int n)
  {
    // Dividing bit flds of the exponent "E" by "m" modulo "PMod":
    return
      (n == 1)
    ? E
    : DivExpRec(E, n, 0UL, 0);
  }

  //=========================================================================//
  // Extraction of Numer and Denom from Zp-Encodings:                        //
  //=========================================================================//
  // Find the Numer and Denom of the minimal total height corresponding to the
  // given Zp representation. Currently this is done just by direct search:
  //
  constexpr int NumerRec(int rep, int height, int denom, int invDenom)
  {
    // Numer is +-(height - denom):
    return
      (((height - denom) * invDenom) % IPMod == rep)
    ? (height - denom)
    : (((IPMod - height + denom) * invDenom) % IPMod == rep)
    ? (denom - height)
    : (denom < height - 1)                                     // Nxt numer != 0
    ? NumerRec(rep, height, denom + 1, InverseModP(denom + 1)) // Rec wrt Denom
    : NumerRec(rep, height + 1, 1, 1);                         // Rec wrt Height
  }

  constexpr int Numer(unsigned long rep)
    // Check for 0, if not, start recursion from height=2, denom=1:
    { return (rep == 0UL) ? 0 : NumerRec(int(rep), 2, 1, 1); }

  constexpr int DenomRec(int rep, int height, int denom, int invDenom)
  {
    // Numer is +-(height - denom):
    return
      ((((height - denom) * invDenom) % IPMod == rep) ||
       (((IPMod - height + denom) * invDenom) % IPMod == rep))
    ? denom
    : (denom < height - 1)                                     // Nxt numer != 0
    ? DenomRec(rep, height, denom + 1, InverseModP(denom + 1)) // Rec wrt Denom
    : DenomRec(rep, height + 1, 1, 1);                         // Rec wrt Height
  }

  constexpr int Denom(unsigned long rep)
    // Check for 0, if not, start recursion from height=2, denom=1:
    { return (rep == 0UL) ? 1 : DenomRec(int(rep), 2, 1, 1); }


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
  // Units corresponding to non-0 Exponent Bit Flds must be the same. This is
  // implemented over templates rather than "constexpr" funcs because we need
  // to provide a clear error message if the units do not unify:
  //
  template<unsigned long e, unsigned long f, unsigned long u, unsigned long v>
  struct UnifyUnitFlds
  {
    static_assert
      (e == 0UL || f == 0UL || u == v, "ERROR: Units do not unify");

    constexpr static unsigned long res =
        (e == 0UL)
      ? v      // "e" is unset, "u" does not matter            -> the res is "v"
      :  u;    // "f" is unset ("v" does not matter) or u == v -> the res is "u"
  };

  template<unsigned long E, unsigned long F, unsigned long U, unsigned long V,
           unsigned long curr, unsigned dim>
  struct UnifyUnitsRec
  {
    constexpr static unsigned long unified =
      UnifyUnitFlds<GetFld(E, dim), GetFld(F, dim),
                    GetFld(U, dim), GetFld(V, dim)>::res;
      
    constexpr static unsigned long res =
      UnifyUnitsRec<E, F, U, V, curr | PutFld(unified, dim), dim + 1>::res;
  };

  template<unsigned long E, unsigned long F, unsigned long U, unsigned long V,
           unsigned long curr>
  struct UnifyUnitsRec<E, F, U, V, curr, NDims>
  {
    constexpr static unsigned long res = curr;
  };

  template<unsigned long E, unsigned long F, unsigned long U, unsigned long V>
  constexpr unsigned long UnifyUnits()
    { return UnifyUnitsRec<E, F, U, V, 0UL, 0>::res; };

  // "UnitsOK": Similar for "UnifyUnits" but for one Exponent Vector. Returns a
  // boolean used in external static assertions, so don't need template recurs-
  // ion here -- can use a recursive constexpr function:
  //
  constexpr bool UnitsOKRec
    (unsigned long E, unsigned long U, unsigned long V, unsigned dim)
  {
    return
      (dim >= NDims)
    ? true                          // All checks finished
    : (GetFld(E,dim) == 0UL || GetFld(U,dim) == GetFld(V,dim))
    ? UnitsOKRec(E, U, V, dim + 1)  // This dim is OK, check next
    : false;                        // This dim failed
  }

  constexpr bool UnitsOK(unsigned long E, unsigned long U, unsigned long V)
    { return UnitsOKRec(E, U, V, 0); }

  //-------------------------------------------------------------------------//
  // Clean-Up (Re-Setting to the default 0) of unused units:                 //
  //-------------------------------------------------------------------------//
  constexpr unsigned long CleanUpFld
    (unsigned long e, unsigned long u, unsigned dim)
  {
    return PutFld((e == 0UL) ? 0UL : u, dim);
  }

  constexpr unsigned long CleanUpUnitsRec
    (unsigned long E, unsigned long U, unsigned long curr, unsigned dim)
  {
    return
      (dim >= NDims)
    ? curr
    : CleanUpUnitsRec
        (E, U, curr | CleanUpFld(GetFld(E,dim), GetFld(U,dim), dim), dim + 1);
  }

  constexpr unsigned long CleanUpUnits(unsigned long E, unsigned long U)
    { return CleanUpUnitsRec(E, U, 0UL, 0); }

  //=========================================================================//
  // Misc:                                                                   //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // Support for Printing Out Exponents and Units:                           //
  //-------------------------------------------------------------------------//
  // FIXME: Make the following functions "constexpr", as they really depend on
  // template params only. For the moment, they require run-time evaluation:
  //
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

  template<int Numer, int Denom>
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
}
// End namespace Bits
}
// End namespace DimTypes
