// vim:ts=2:et
//===========================================================================//
//                          "DimTypes/Bits/Macros.h":                        //
//                   Macros for Declaring Dimensions and Units               //
//===========================================================================//
#pragma  once
#include <cstdio>
#include <cassert>

//===========================================================================//
// MACROS for User-Level Declaration of Dimensions and Units:                //
//===========================================================================//
// They can be invoked in any namespace where the corresp declarations need to
// be placed. All names created by them are placed in the macro expansion con-
// text, NOT in the "DimTypes" namespace:
//===========================================================================//
// "FOR_EACH": Recursive Macro Application:                                  //
//===========================================================================//
// From: https://www.scs.stanford.edu/~dm/blog/va-opt.html
// Recursively applies macro "ACTION" to the arg list, up to ~30 entries long:
//
#ifdef  PARENS
#undef  PARENS
#endif
#ifdef  EXPAND
#undef  EXPAND
#endif
#ifdef  EXPAND1
#undef  EXPAND1
#endif
#ifdef  EXPAND2
#undef  EXPAND2
#endif
#ifdef  EXPAND3
#undef  EXPAND3
#endif
#ifdef  FOR_EACH_HELPER
#undef  FOR_EACH_HELPER
#endif
#ifdef  FOR_EACH_AGAIN
#undef  FOR_EACH_AGAIN
#endif
#ifdef  FOR_EACH
#undef  FO_EACH
#endif

#define PARENS ()

#define EXPAND( ...) EXPAND1(EXPAND1(EXPAND1(__VA_ARGS__)))
#define EXPAND1(...) EXPAND2(EXPAND2(EXPAND2(__VA_ARGS__)))
#define EXPAND2(...) EXPAND3(EXPAND3(EXPAND3(__VA_ARGS__)))
#define EXPAND3(...) __VA_ARGS__

#define FOR_EACH(action, ...)                                    \
  __VA_OPT__(EXPAND(FOR_EACH_HELPER(ACTION, __VA_ARGS__)))
#define FOR_EACH_HELPER(action, arg1, ...)                       \
  action(arg1)                                                   \
  __VA_OPT__(FOR_EACH_AGAIN PARENS (action, __VA_ARGS__))
#define FOR_EACH_AGAIN() FOR_EACH_HELPER

//===========================================================================//
// "DECLARE_DIMS":                                                           //
//===========================================================================//
// Format:
// DECLARE_DIMS(
//   RepT,
//   FundUnitName
//   [,(AnotherUnit, AnotherUnitValInFundUnits)...]
// )
#ifdef  DECLARE_DIMS
#undef  DECLARE_DIMS
#endif
#define DECLARE_DIMS(...) \
  /*-----------------------------------------------------------------------*/ \
  /* "DimsE": Enum of all Dims to be used.                                 */ \
  /* Their encodings are assigned automatically, starting from 0:          */ \
  /*-----------------------------------------------------------------------*/ \
  enum class DimsE: unsigned { __VA_ARGS__ }; \
  \
  /*-----------------------------------------------------------------------*/ \
  /* Template prototypes for subsequent specialisation                     */ \
  /* (with Dims and Units encodings):                                      */ \
  /*-----------------------------------------------------------------------*/ \
  template<unsigned Dim, unsigned Unit, typename RepT> \
  constexpr RepT UnitScale        = RepT(1.0); \
  \
  template<unsigned Dim, unsigned Unit>      \
  constexpr char const* UnitName  = nullptr; \
  \
  /*-----------------------------------------------------------------------*/ \
  /* Output Function for "DimQs": NB: it is NOT "constexpr"!               */ \
  /*-----------------------------------------------------------------------*/ \
  template<unsigned long E, unsigned long U, typename RepT> \
  char* Put \
  ( \
    DimTypes::DimQ<E, U, RepT> a_dimq, \
    char*                      a_buff, \
    char const*                a_end   \
  ) \
  { \
    assert(a_buff != nullptr && a_buff < a_end); \
    char* curr = a_buff;               \
    int   N    = int(a_end - a_buff);  \
    \
    /* First, the Magnitude. As "RepT" can be any type, we apply some rules */ \
    /* of thumb here,  trying to fugure out whether it is a real or complex */ \
    /* one, and to choose the formats accordingly:                          */ \
    RepT  mag  = a_dimq.Magnitude();   \
    if constexpr(DimTypes::Bits::IsComplex<RepT>) \
    { \
      double  magRe  = double(mag.real()); \
      double  magIm  = double(mag.imag()); \
      \
      curr += snprintf(buff,  N, "(%.16e %c %.16e * I)", \
                       magRe, (magIm < 0.0) ? '-' : '+', fabs(magIm)); \
    } \
    else \
      curr += snprintf(buff,  N, "%.16e", double(mag));  \
    \
    /* Now the Units and Exponents: */ \
    for (unsigned dim = 0; dim < DimTypes::NDims; ++dim) \
    { \
      /* NB: the unit name is returned as a string by UnitName<Dim,Unit>    */ \
      /* which is specialised for the types actually used:                  */ \
      char const* unitName = UnitName<Dim, DimTypes::Bits::GetFld(U,Dim)>;  \

    constexpr auto NumDen  = \
      DimTypes::Bits::GetNumerAndDenom (DimTypes::Bits::GetFld(E,Dim));   \
    constexpr int     Numer  = NumDen.first;  \
    constexpr unsigned Denom = NumDen.second; \
    curr += DimTypes::Bits::UnitStr<Numer, Denom>::put

            (curr, UnitImage<Dim,DimTypes::Bits::GetFld(U,Dim)>);  \
    PUT_UNIT_STR(0) \
    PUT_UNIT_STR(1) \
    PUT_UNIT_STR(2) \
    PUT_UNIT_STR(3) \
    PUT_UNIT_STR(4) \
    PUT_UNIT_STR(5) \
    PUT_UNIT_STR(6) \
  } \
  \



//===========================================================================//
// "DECLARE_DIM_UNITS":                                                      //
//===========================================================================//
// Declaring all Units for a given "DimName" (which is a member of "DimsE") andr
// thus assigned per-dimension encodings to them.
// Automatically invokes "DECLARE_UNIT" for the 1st (Fundamental) unit, but for
// all other (aux) units listed here, "DECLARE_UNIT" must subsequently be invok-
// ed explicitly:
//
#ifdef  DECLARE_DIM_UNITS
#undef  DECLARE_DIM_UNITS
#endif
#define DECLARE_DIM_UNITS(DimName, RepT, FundUnitName, ...)       \
  static_assert(unsigned(DimsE::DimName) < DimTypes::Bits::NDims, \
                "Too many Dimensions"); \
  enum class DimName##UnitsE: int  { FundUnitName, __VA_ARGS__ }; \
  DECLARE_UNIT(DimName, RepT, FundUnitName, 1.0)

//===========================================================================//
// "DECLARE_UNIT":                                                           //
//===========================================================================//
// For a given unit, declares a "DimQ" short-cut, a stringifier specialisation,
// the scale (value in Fundamental Units) specialn, and a conversion funcs for
// arbitrary "DimQ"s into this Unit.
// NB:
// (*) With the naming conventioned used, it is OK to have same-named units for
//     different dims;
// (*) Modular arithmetic is not used for Units, so any Unit codes up to and
//     including "PMask" are OK:
//
#ifdef  DECLARE_UNIT
#undef  DECLARE_UNIT
#endif
#define DECLARE_UNIT(DimName, RepT, UnitName, UnitVal)  \
  static_assert(unsigned(DimName##UnitsE::UnitName) <=  \
                DimTypes::Bits::PMask, "Too many Units for this Dim"); \
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
  /* NB: It has the value of 1.0 in ITSELF; NOT converted automatically to */ \
  /*     the corresp Fundamental Units:                                    */ \
  constexpr DimName##_##UnitName##_##T DimName##_##UnitName(RepT(1.0)); \
  \
  /*-----------------------------------------------------------------------*/ \
  /* UnitScale specialisation for this Unit:                               */ \
  /*-----------------------------------------------------------------------*/ \
  /* HERE the Unit Conversion Factor (UnitSCale) is actually introduced:   */ \
  template<> \
  constexpr RepT UnitScale  \
      <unsigned(DimsE::DimName), \
       unsigned(DimName##UnitsE::UnitName), RepT> = RepT(UnitVal); \
  \
  /*-----------------------------------------------------------------------*/ \
  /* UnitName specialisation for this Unit:                                */ \
  /*-----------------------------------------------------------------------*/ \
  template<> \
  char const* UnitImage \
    <unsigned(DimsE::DimName),            \
     unsigned(DimName##UnitsE::UnitName)> = #UnitName;  \
  \
  /*-----------------------------------------------------------------------*/ \
  /* Conversion function for DimQs: To_DimName_UnitName:                   */ \
  /*-----------------------------------------------------------------------*/ \
  template<unsigned long E, unsigned long U>      \
  constexpr DimTypes::DimQ   \
    <E, \
     DimTypes::Bits::SetUnit \
       (U, unsigned(DimsE::DimName), unsigned(DimName##UnitsE::UnitName)), \
     RepT> \
  To_##DimName##_##UnitName(DimTypes::DimQ<E,U,RepT> a_dimq) \
  { \
    constexpr unsigned Dim       = unsigned(DimsE::DimName); \
    constexpr unsigned OldUnit   = DimTypes::Bits::GetFld     (U, Dim);    \
    constexpr unsigned NewUnit   = unsigned(DimName##UnitsE::UnitName);    \
    constexpr auto     ExpNumDen = \
      DimTypes::Bits::GetNumerAndDenom(DimTypes::Bits::GetFld(E, Dim));    \
    constexpr int      Numer     = ExpNumDen.first;  \
    constexpr unsigned Denom     = ExpNumDen.second; \
    return \
      DimTypes::DimQ<E, DimTypes::Bits::SetUnit(U, Dim, NewUnit),  RepT>   \
      ( \
        a_dimq.Magnitude() * \
        DimTypes::Bits::FracPower<Numer, Denom, RepT> \
          /* OldScale / NewScale: */                  \
          (UnitScale<Dim, OldUnit, RepT> / UnitScale<Dim, NewUnit, RepT>)  \
      ); \
  }

