// vim:ts=2:et
//===========================================================================//
//                          "DimTypes/DimTypes.hpp":                         //
//          Types for Physical Dimensions and Dimensioned Quantities         //
//===========================================================================//
#pragma  once

#include "Bits/FracPower.hpp"
#include "Bits/Encodings.hpp"
#include <complex>
#include <cassert>

namespace DimTypes
{
  using namespace Bits;

  //=========================================================================//
  // "DimQ": The External Dimensioned Type:                                  //
  //=========================================================================//
  // "DimQ" is actually a Field "T" statically parameterised by the encoded
  // Dimension Exponent "E" and the vector of units "U":
  // NB: "RepT" and "DimQ" are everywhere passed BY COPY (except for result vals
  // of in-place update operators),  which is probably optimal even if "RepT" is
  // a complex or multiple-precision type, to provide better data locality:
  //
  template<unsigned long E, unsigned long U, typename RepT = double>
  class DimQ
  {
  private:
    RepT m_val;  // The actual value (magnitude)
  
  public:
    //=======================================================================//
    // Ctors and Accessors:                                                  //
    //=======================================================================//
    // Non-Default Ctor: Lifting a value into "DimQ" (by default, 0):
    //
    constexpr explicit DimQ(RepT val = RepT(0.0)): m_val(val) {}
  
    // Copy Ctor:
    constexpr DimQ(DimQ const& a_right): m_val(a_right.m_val) {}

    // Assignment:
    constexpr DimQ& operator= (DimQ const& a_right)
    {
      m_val = a_right.m_val;
      return *this;
    }

    // The dimensioned unit for the given dimensioned quantity. NB: This method
    // is not a ctor: it takes an existing "DimQ" and returns a similar one but
    // with the unitary magnitude:
    //
    constexpr DimQ UnitOf() const  { return DimQ(RepT(1.0));  }
  
    // Getting the exponent code (this is primarily for testing):
    constexpr unsigned long        GetDimsCode ()    const { return E; }
  
    // Getting the Units code (also for testing):
    constexpr unsigned long        GetUnitsCode()    const { return U; }
  
    // The Magnitude of the dimensioned quantity (ie its value expressed in the
    // corresp dimensioned units). The magnitude is by itself dimension-less.
    // XXX: direct access to the dimension-less magnitude  may seem  to be an
    // unsafe feature but it could be emulated anyway (by dividing the qty by
    // its unit), so we provide direct access to it for efficiency:
    //
    constexpr RepT Magnitude() const  { return m_val; }
  
    // The following function converts "DimQ" to the Underlying Field Type
    // ("RepT"), but ONLY if the arg is actually dimension-less:
    constexpr explicit operator RepT() const
    {
      static_assert(E == 0, "Must be a DimLess qty");
      return m_val;
    }
  
    //=======================================================================//
    // Arithmetic Operations -- Dimensions Unchanged:                        //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Addition and Subtraction of "DimQ"s:                                  //
    //-----------------------------------------------------------------------//
    // This is only possible if the dimensions and rep are same, and the units
    // can be unified (i.e. same units for non-0 exponents), in which case the
    // 1st arg's units can still be used:
    //
    // Addition:
    template<unsigned long V>
    constexpr DimQ operator+   (DimQ<E,V> a_right) const
    {
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify");
      return DimQ(m_val + a_right.Magnitude());
    }
  
    template<unsigned long V>
    constexpr DimQ& operator+= (DimQ<E,V> a_right)
    {
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify");
      m_val += a_right.Magnitude();
      return *this;
    }
  
    // Subtraction of "DimQs":
    // Same constraints as for addition:
    template<unsigned long V>
    constexpr DimQ operator- (DimQ<E,V,RepT> a_right) const
    {
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify");
      return DimQ(m_val - a_right.Magnitude());
    }
  
    template<unsigned long V>
    constexpr DimQ& operator-=(DimQ<E,V,RepT> a_right)
    {
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify");
      m_val -= a_right.Magnitude();
      return *this;
    }
  
    //-----------------------------------------------------------------------//
    // Multiplication and Division by a Dimension-Less Factor:               //
    //-----------------------------------------------------------------------//
    // NB: the following method can be used to create new dimensioned values
    // from the fundamental units. XXX: when used this way, it is slighly in-
    // efficient because of multiplication by 1 and extra copying;   but the
    // impact of that should be negligible:
    //
    // Multiplication:
    constexpr DimQ operator*(RepT a_right) const
      { return DimQ(m_val * a_right); }
  
    constexpr DimQ& operator*= (RepT a_right)
    {
      m_val *= a_right;
      return *this;
    }
  
    constexpr friend DimQ operator* (RepT a_left, DimQ a_right)
      { return DimQ(a_left * a_right.Magnitude()); }
  
    // In particular, Unary Negation:
    constexpr  DimQ operator- () const
      { return DimQ(- m_val); }
  
    // Division by a Dimension-Less Factor:
    constexpr  DimQ operator/   (RepT a_right) const
      { return DimQ(m_val / a_right); }
  
    constexpr  DimQ& operator/= (RepT a_right)
    {
      m_val /= a_right;
      return *this;
    }
  
    //-----------------------------------------------------------------------//
    // "Abs"; "Floor", "Ceil", "Round" (if supported by "RepT"):             //
    //-----------------------------------------------------------------------//
    constexpr DimQ Abs  () const { return DimQ(Bits::Abs  (m_val)); }
    constexpr DimQ Floor() const { return DimQ(Bits::Floor(m_val)); }
    constexpr DimQ Ceil () const { return DimQ(Bits::Ceil (m_val)); }
    constexpr DimQ Round() const { return DimQ(Bits::Round(m_val)); }

    //=======================================================================//
    // Arithmetic operations which result in Dimensions change:              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Multiplication and Division of "DimQ"s:                               //
    //-----------------------------------------------------------------------//
    // Dimension exponents are added (or subtracted) up; units must unify;
    // afterwards, units of zero-exp dims are reset, otherwise false type-che-
    // cking failures may occur:
    //
    template<unsigned long F, unsigned long V>
    constexpr DimQ<AddExp(E,F),
                   CleanUpUnits(AddExp(E,F), UnifyUnits(E,F,U,V)), RepT>
    operator* (DimQ<F,V,RepT> a_right) const
    {
      return DimQ<AddExp(E,F),
                  CleanUpUnits(AddExp(E,F),  UnifyUnits(E,F,U,V)), RepT>
             (m_val * a_right.Magnitude());
    }
  
    // Division of "DimQ"s:
    // Dimension exponents are subtracted. Same treatment of units as for mult:
    template<unsigned long F, unsigned long V>
    constexpr DimQ<SubExp(E,F),
                   CleanUpUnits(SubExp(E,F), UnifyUnits(E,F,U,V)), RepT>
    operator/ (DimQ<F,V,RepT> a_right) const
    {
      return DimQ<SubExp(E,F),
                  CleanUpUnits(SubExp(E,F),  UnifyUnits(E,F,U,V)), RepT>
             (m_val / a_right.Magnitude());
    }
  
    constexpr friend DimQ<SubExp(0UL,E),U,RepT>
    operator/ (RepT a_left, DimQ a_right)
      { return DimQ<SubExp(0UL,E),U,RepT>(a_left / a_right.Magnitude()); }
  
    // Integral and Rational Powers -- the power must be given by constant exprs
    // so it is analysable at compile time:
    // NB: for applying these methods from outside "DimQ",  explicit  "template
    // pow" expr is required, otherwise the template is not syntactically recog-
    // nised!
    // NB: "CleanUpUnits" is required in case if M==0:
    //
    // "IPow": Integral Power:
    template<int M>
    constexpr DimQ<MultExp(E,M), CleanUpUnits(MultExp(E,M), U), RepT>
    IPow() const
    {
      return
        DimQ<MultExp(E,M), CleanUpUnits(MultExp(E,M), U), RepT>
            (IntPower<RepT, M>::res(m_val));
    }

    // "RPow": (General) Rational Power:
    template<int M, int N>
    DimQ<DivExp(MultExp(E,M),N),
         CleanUpUnits(DivExp(MultExp(E,M),N), U), RepT>
    RPow() const
    {
      // Multiples of "PMod" (including 0) are uninvertible in the curr rep;
      // and for a fractional power, "constexpr" is not possible yet:
      //
      static_assert(N > 0 && N % IPMod != 0, "RPow: Invalid Denom");
      return
        DimQ<DivExp(MultExp(E,M),N),
             CleanUpUnits(DivExp(MultExp(E,M),N), U), RepT>
            (FracPower<M,N,RepT>(m_val));
    }

    // Shortcuts: "SqRt" and "CbRt". Here M==1, so "CleanUpUnits" is not
    // required:
    constexpr DimQ<DivExp(E,2),U,RepT> SqRt() const { return RPow<1,2>(); }
    constexpr DimQ<DivExp(E,3),U,RepT> CbRt() const { return RPow<1,3>(); }

    //-----------------------------------------------------------------------//
    // Comparison operators:                                                 //
    //-----------------------------------------------------------------------//
    // Require same representation, same dimensions and unifyable units.
    // NB:
    // (*) Inequalities would result in a compile-time error if not implemented
    //     for the corresp underlying type (eg a Complex) but only at the use
    //     point!
    // (*) Interestingly, a space is allowed between e.g. "operator" and "==" :
    //
#   ifdef  DIMQ_CMP
#   undef  DIMQ_CMP
#   endif
#   define DIMQ_CMP(Op) \
    template<unsigned long V> \
    constexpr bool operator Op (DimQ<E,V,RepT> a_right) const \
    { \
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify"); \
      return m_val Op a_right.Magnitude(); \
    }
    DIMQ_CMP(==)
    DIMQ_CMP(!=)
    DIMQ_CMP(>)
    DIMQ_CMP(>=)
    DIMQ_CMP(<)
    DIMQ_CMP(<=)
  
    constexpr bool IsZero  () const  { return m_val == RepT(0.0); }
    constexpr bool IsFinite() const  { return IsFinite(m_val); }

    // NB: The following methods will not compile if the field "RepT" is not
    // ordered (eg for complex numbers):
    constexpr bool IsNeg   () const  { return m_val <  RepT(0.0); }
    constexpr bool IsPos   () const  { return m_val >  RepT(0.0); }
  };

  //=========================================================================//
  // Syntactic Sugar: Above methods in Prefix notation:                      //
  //=========================================================================//
  // XXX: INCREDIBLY, in C++20, these functions are visible directly (ie w/o the
  // "Dims::" prefix!!!) outside the "Dims" namespace, if the corresp arg "DimQ"
  // obj is there. So their visibility is very much like that of "DimQ" methods:
  //
  template<unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<E,U,RepT> UnitOf(DimQ<E,U,RepT> a_dimq)
    { return a_dimq.UnitOf(); }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr unsigned long GetDimsCode (DimQ<E,U,RepT> a_dimq)
    { return a_dimq.GetDimsCode();  }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr unsigned long GetUnitsCode(DimQ<E,U,RepT> a_dimq)
    { return a_dimq.GetUnitsCode(); }
  
  template<unsigned long E, unsigned long U, typename RepT>
  constexpr RepT Magnitude(DimQ<E,U,RepT> a_dimq)
    { return a_dimq.Magnitude();    }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<E,U,RepT>  Abs  (DimQ<E,U,RepT> a_right)
    { return a_right.Abs();    }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<E,U,RepT>  Floor(DimQ<E,U,RepT> a_right)
    { return a_right.Floor();  }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<E,U,RepT>  Ceil (DimQ<E,U,RepT> a_right)
    { return a_right.Ceil();  }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<E,U,RepT>  Round(DimQ<E,U,RepT> a_right)
    { return a_right.Round(); }

  template<int M, unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<MultExp(E,M), CleanUpUnits(MultExp(E,M),U), RepT>
  IPow(DimQ<E,U,RepT> a_right)
    { return a_right.template IPow<M>(); }

  template<int M, int N, unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<DivExp(MultExp(E,M),N),
                 CleanUpUnits(DivExp(MultExp(E,M),N), U), RepT>
  RPow(DimQ<E,U,RepT> a_right)
    { return a_right.template RPow<M,N>(); }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<DivExp(E,2),U,RepT> SqRt(DimQ<E,U,RepT> a_right)
    { return a_right.SqRt(); }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr DimQ<DivExp(E,3),U,RepT> CbRt(DimQ<E,U,RepT> a_right)
    { return a_right.CbRt(); }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr bool IsZero  (DimQ<E,U,RepT> a_right)
    { return a_right.IsZero();     }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr bool IsFinite(DimQ<E,U,RepT> a_right)
    { return a_right.IsFinite();   }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr bool IsNeg   (DimQ<E,U,RepT> a_right)
    { return a_right.IsNegative(); }

  template<unsigned long E, unsigned long U, typename RepT>
  constexpr bool IsPos   (DimQ<E,U,RepT> a_right)
    { return a_right.IsPositive(); }
}
// End namespace DimTypes

//===========================================================================//
// Macros for User-Level Declarations of Dimensions and Units:               //
//===========================================================================//
// They can be invoked in any namespace where the corresp declarations need to
// be placed:
//---------------------------------------------------------------------------//
// "DECLARE_DIMS":                                                           //
//---------------------------------------------------------------------------//
// First macro to be invoked at the user level. Declares the enum type (called
// "DimsE") of all dimensions to be used. This assigns encodings to dimensions,
// starting from 0:
//
#ifdef  DECLARE_DIMS
#undef  DECLARE_DIMS
#endif
#define DECLARE_DIMS(...) \
  enum class DimsE: unsigned { __VA_ARGS__ }; \
  \
  template<unsigned Dim, unsigned Unit, typename RepT> \
  RepT UnitScale; \
  \
  template<unsigned Dim, unsigned Unit> \
  char const* UnitImage;

//---------------------------------------------------------------------------//
// "DECLARE_DIM_UNITS":                                                      //
//---------------------------------------------------------------------------//
// Declaring all Units of a given Dimension and thus assigned per-dimension en-
// codings to them. Automatically invokes "DECLARE_UNIT" for the 1st (Fundamen-
// tal) unit, but for all other (aux) units listed here,   "DECLARE_UNIT" must
// subsequently be invoked explicitly:
//
#ifdef  DECLARE_DIM_UNITS
#undef  DECLARE_DIM_UNITS
#endif
#define DECLARE_DIM_UNITS(DimName, RepT, FundUnitName, ...)       \
  static_assert(unsigned(DimsE::DimName) < DimTypes::Bits::NDims, \
                "Too many Dimensions"); \
  enum class DimName##UnitsE: int  { FundUnitName, __VA_ARGS__ }; \
  DECLARE_UNIT(DimName, RepT, FundUnitName, 1.0)

//---------------------------------------------------------------------------//
// "DECLARE_UNIT":                                                           //
//---------------------------------------------------------------------------//
// For a given unit, declares a factory function,   the stringifier specialsn,
// the scale (value in Fundamental Units) specialn, and a conversion funcs for
// arbitrary "DimQ"s into this Unit.
// NB: with the naming conventioned used, it is OK to have same-named units for
// different dims:
//
#ifdef  DECLARE_UNIT
#undef  DECLARE_UNIT
#endif
#define DECLARE_UNIT(DimName, RepT, UnitName, UnitVal)    \
  static_assert(unsigned(DimName##UnitsE::UnitName) <= \
                DimTypes::Bits::PMod, "Too many Units for this Dim"); \
  \
  /*-----------------------------------------------------------------------*/ \
  /* DimName_UnitName_T as a type:                                         */ \
  /*-----------------------------------------------------------------------*/ \
  using DimName##_##UnitName##_##T = \
    DimTypes::DimQ \
      <DimTypes::Bits::DimExp(unsigned(DimsE::DimName)), \
       DimTypes::Bits::MkUnit(unsigned(DimsE::DimName),  \
                              unsigned(DimName##UnitsE::UnitName)), \
       RepT>; \
  \
  /*-----------------------------------------------------------------------*/ \
  /* DimName_UnitName, as a templated const:                               */ \
  /*-----------------------------------------------------------------------*/ \
  /* NB: It has the value of 1 in ITSELF; NOT converted automatically to   */ \
  /*     the corresp Fundamental Units:                                    */ \
  constexpr DimName##_##UnitName##_##T DimName##_##UnitName(RepT(1.0)); \
  \
  /*-----------------------------------------------------------------------*/ \
  /* UnitScale specialisation for this Unit:                               */ \
  /*-----------------------------------------------------------------------*/ \
  /* HERE the Conversion Factor is actually introduced:                    */ \
  template<> \
  constexpr RepT UnitScale  \
      <unsigned(DimsE::DimName), \
       unsigned(DimName##UnitsE::UnitName), RepT> = RepT(UnitVal); \
  \
  /*-----------------------------------------------------------------------*/ \
  /* UnitImage specialisation for this Unit:                               */ \
  /*-----------------------------------------------------------------------*/ \
  template<> \
  char const* UnitImage \
    <unsigned(DimsE::DimName),            \
     unsigned(DimName##UnitsE::UnitName)> = #UnitName;  \
  \
  /* Conversion function for DimQs: To_DimName_UnitName: */    \
  template<unsigned long E, unsigned long U>      \
  constexpr DimTypes::DimQ   \
    <E, \
     DimTypes::Bits::SetUnit \
       (U, unsigned(DimsE::DimName), unsigned(DimName##UnitsE::UnitName)), \
     RepT> \
  To_##DimName##_##UnitName(DimTypes::DimQ<E,U,RepT> a_dimq) \
  { \
    constexpr unsigned Dim     = unsigned(DimsE::DimName); \
    constexpr unsigned OldUnit = DimTypes::Bits::GetFld     (U, Dim);     \
    constexpr unsigned NewUnit = unsigned(DimName##UnitsE::UnitName);     \
    constexpr auto     NumDen  = \
      DimTypes::Bits::GetNumerAndDenom(DimTypes::Bits::GetFld(E, Dim));   \
    constexpr int      Numer   = NumDen.first;  \
    constexpr unsigned Denom   = NumDen.second; \
    return \
      DimTypes::DimQ<E, DimTypes::Bits::SetUnit(U, Dim, NewUnit),  RepT>  \
      ( \
        a_dimq.Magnitude() * \
        DimTypes::Bits::FracPower<Numer, Denom, RepT> \
          /* OldScale / NewScale: */                  \
          (UnitScale<Dim, OldUnit, RepT> / UnitScale<Dim, NewUnit, RepT>) \
      ); \
  }

//---------------------------------------------------------------------------//
// "DECLARE_DIM_STR":                                                        //
//---------------------------------------------------------------------------//
// A macro which generates (in the context where it is invoked, and where the
// Dims have been declared) the "ToStr" function for "DimQ"s.  XXX: This func
// is NOT thread-safe and is quite inefficient, so it should only be used for
// the final output or for debugging:
//
#ifdef  PUT_UNIT_STR
#undef  PUT_UNIT_STR
#endif
#define PUT_UNIT_STR(Dim) \
  { \
    constexpr auto NumDen  = \
      DimTypes::Bits::GetNumerAndDenom (DimTypes::Bits::GetFld(E,Dim));   \
    constexpr int     Numer  = NumDen.first;  \
    constexpr unsigned Denom = NumDen.second; \
    curr += DimTypes::Bits::UnitStr<Numer, Denom>::put \
            (curr, UnitImage<Dim,DimTypes::Bits::GetFld(U,Dim)>);  \
  }

#ifdef  DECLARE_DIM_STR
#undef  DECLARE_DIM_STR
#endif
#define DECLARE_DIM_STR   \
  template<typename T> \
  constexpr bool IsComplex = false; \
  template<typename T>    \
  constexpr bool IsComplex<std::complex<T>> = true; \
  \
  template<unsigned long E, unsigned long U, typename RepT> \
  char const* ToStr(DimTypes::DimQ<E, U, RepT> a_dimq) \
  { \
    constexpr int  N = 256;  \
    static char    buff[N];  \
    char*  curr =  buff;     \
    /* First, the Magnitude. As "RepT" can be any type, we apply some rules */ \
    /* of thumb here,  trying to fugure out whether it is a real or complex */ \
    /* one, and to choose the formats accordingly:                          */ \
    RepT mag = a_dimq.Magnitude(); \
    if constexpr(IsComplex<RepT>)  \
    { \
      double  magRe  = double(mag.real()); \
      double  magIm  = double(mag.imag()); \
      curr += snprintf(buff,  N, "(%.16e %c %.16e * I)", \
                       magRe, (magIm < 0.0) ? '-' : '+', fabs(magIm)); \
    } \
    else \
      curr += snprintf(buff,  N, "%.16e", double(mag)); \
    \
    /* NB: the unit name is returned by UnitImage<Dim,Unit> which is speci- */ \
    /* alised for the types actually used:                                  */ \
    PUT_UNIT_STR(0) \
    PUT_UNIT_STR(1) \
    PUT_UNIT_STR(2) \
    PUT_UNIT_STR(3) \
    PUT_UNIT_STR(4) \
    PUT_UNIT_STR(5) \
    PUT_UNIT_STR(6) \
    assert(curr < buff + N && *curr == '\0'); \
    return buff;    \
  }
// EOF
