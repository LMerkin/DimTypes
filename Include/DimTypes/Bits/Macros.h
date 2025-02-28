// vim:ts=2:et
//===========================================================================//
//                          "DimTypes/Bits/Macros.h":                        //
//                   Macros for Declaring Dimensions and Units               //
//===========================================================================//
#pragma  once
#include <cstdio>
#include <cstring>
#include <cassert>
#include <stdexcept>
#include <array>
#include <format>

//===========================================================================//
// MACROS for User-Level Declaration of Dimensions and Units:                //
//===========================================================================//
// They can be invoked in any namespace where the corresp declarations need to
// be placed. All names created by them are placed in the macro expansion con-
// text, NOT in the "DimTypes" namespace:
//
//===========================================================================//
// "FOR_EACH_DIM":                                                           //
//===========================================================================//
// From: https://www.scs.stanford.edu/~dm/blog/va-opt.html
// Iteratively applies macro "Action" to the Args list, up to ~30 entries long:
//
#ifdef  PARENS
#undef  PARENS
#endif

#ifdef  EXPAND_DIMS
#undef  EXPAND_DIMS
#endif

#ifdef  EXPAND_DIMS1
#undef  EXPAND_DIMS1
#endif

#ifdef  EXPAND_DIMS2
#undef  EXPAND_DIMS2
#endif

#ifdef  EXPAND_DIMS3
#undef  EXPAND_DIMS3
#endif

#ifdef  FOR_EACH_DIM_AGAIN
#undef  FOR_EACH_DIM_AGAIN
#endif

#ifdef  FOR_EACH_DIM_HELPER
#undef  FOR_EACH_DIM_HELPER
#endif

#ifdef  FOR_EACH_DIM
#undef  FOR_EACH_DIM
#endif

#define PARENS ()

#define EXPAND_DIMS(...) \
        EXPAND_DIMS1(EXPAND_DIMS1(EXPAND_DIMS1(__VA_ARGS__)))

#define EXPAND_DIMS1(...) \
        EXPAND_DIMS2(EXPAND_DIMS2(EXPAND_DIMS2(__VA_ARGS__)))

#define EXPAND_DIMS2(...) \
        EXPAND_DIMS3(EXPAND_DIMS3(EXPAND_DIMS3(__VA_ARGS__)))

#define EXPAND_DIMS3(...) __VA_ARGS__

#define FOR_EACH_DIM(Action, ...)                   \
  __VA_OPT__(EXPAND_DIMS(FOR_EACH_DIM_HELPER(Action, __VA_ARGS__)))

#define FOR_EACH_DIM_HELPER(Action, Arg1, ...)      \
  /* NB: Arg1 is assumed to be enclosed in ()s! */  \
  Action(Arg1)                                      \
  __VA_OPT__(FOR_EACH_DIM_AGAIN PARENS (Action, __VA_ARGS__))

#define FOR_EACH_DIM_AGAIN() FOR_EACH_DIM_HELPER

//===========================================================================//
// "FOR_EACH_OTHER_UNIT":                                                    //
//===========================================================================//
// Iterative "Action" application, similar to "FOR_EACH_DIM". Named separately
// in order to avoid errors when the former macro is invoked  recursively from
// the latter):
//
#ifdef  EXPAND_UNITS
#undef  EXPAND_UNITS
#endif

#ifdef  EXPAND_UNITS1
#undef  EXPAND_UNITS1
#endif

#ifdef  EXPAND_UNITS2
#undef  EXPAND_UNITS2
#endif

#ifdef  EXPAND_UNITS3
#undef  EXPAND_UNITS3
#endif

#ifdef  FOR_EACH_OTHER_UNIT_AGAIN
#undef  FOR_EACH_OTHER_UNIT_AGAIN
#endif

#ifdef  FOR_EACH_OTHER_UNIT_HELPER
#undef  FOR_EACH_OTHER_UNIT_HELPER
#endif

#ifdef  FOR_EACH_OTHER_UNIT
#undef  FOR_EACH_OTHER_UNIT
#endif

#define EXPAND_UNITS(...) \
        EXPAND_UNITS1(EXPAND_UNITS1(EXPAND_UNITS1(__VA_ARGS__)))

#define EXPAND_UNITS1(...) \
        EXPAND_UNITS2(EXPAND_UNITS2(EXPAND_UNITS2(__VA_ARGS__)))

#define EXPAND_UNITS2(...) \
        EXPAND_UNITS3(EXPAND_UNITS3(EXPAND_UNITS3(__VA_ARGS__)))

#define EXPAND_UNITS3(...) __VA_ARGS__

#define FOR_EACH_OTHER_UNIT_AGAIN() FOR_EACH_OTHER_UNIT_HELPER

#define FOR_EACH_OTHER_UNIT(Action, DimName, ...)   \
  __VA_OPT__(EXPAND_UNITS(FOR_EACH_OTHER_UNIT_HELPER(Action, DimName, \
                                                     __VA_ARGS__)))

#define FOR_EACH_OTHER_UNIT_HELPER(Action, DimName, Arg1, ...)  \
  /* NB: Arg1 is assumed to be enclosed in ()s! */  \
  Action(DimName, Arg1)                         \
  __VA_OPT__(FOR_EACH_OTHER_UNIT_AGAIN PARENS   \
            (Action,  DimName,  __VA_ARGS__))

//===========================================================================//
// "DECLARE_DIMS":                                                           //
//===========================================================================//
// Usage:
// DECLARE_DIMS(
//   RepT, [MaxDims],
//   (FundUnitName [,(AnotherUnitName, AnotherUnitValInFundUnits)...]),
//   ...
// )
// For a given Unit, declares a "DimQ" short-cut type, a stringifier speciali-
// sation, the scale (value in the corresp Fundamental Unit) specialisation,
// and a conversion function for arbitrary "DimQ"s into this Unit.
// NB:
// (*) "MaxDims" param may be empty, but the comma following it is mandatory;
//     if specified explicitly, "MaxDims" may be 7, 8, or 9;
// (*) with the naming conventions used, it is OK to have same-named Units for
//     different Dims, though this is very rarely needed;
// (*) modular arithmetic is not used for Units,  so any Unit codes up to and
//     including the corresp "PMask" are OK:
//
#ifdef  DECLARE_DIMS
#undef  DECLARE_DIMS
#endif
#define DECLARE_DIMS(RepT, MaxDims, ...) \
  /*-----------------------------------------------------------------------*/ \
  /* "DimQ_MaxDims", "DimQ_RepT" and "DimQ_Encs":                          */ \
  /*-----------------------------------------------------------------------*/ \
  constexpr static unsigned   DimQ_MaxDims = \
    ( unsigned{MaxDims} == 0 )  \
    ? DimTypes::Bits::DefMaxDims \
    : unsigned{MaxDims};   \
  using                       DimQ_RepT    = RepT; \
  using                       DimQ_Encs    = \
    DimTypes::Bits::Encodings<RepT, DimQ_MaxDims>; \
  \
  /*-----------------------------------------------------------------------*/ \
  /* "DimsE": Enum of all Dims to be used.                                 */ \
  /* Their encodings are assigned automatically, starting from 0:          */ \
  /*-----------------------------------------------------------------------*/ \
  enum class DimsE: unsigned \
  {  \
    FOR_EACH_DIM(GET_DIM_NAME_COMMA, __VA_ARGS__) \
  }; \
  \
  /*-----------------------------------------------------------------------*/ \
  /* "UnitNameStr" and "UnitScale" template prototypes for subsequent      */ \
  /* specialisation with actual Dims and Units encodings:                  */ \
  /*-----------------------------------------------------------------------*/ \
  template<unsigned Dim, unsigned Unit> \
  char const* UnitNameStr   = nullptr;  \
  \
  template<unsigned Dim, unsigned Unit> \
  DimQ_RepT UnitScale = DimQ_RepT(NAN); \
  \
  /*-----------------------------------------------------------------------*/ \
  /* For each Dim, declare all its Units, their vals and conversion funcs: */ \
  /*-----------------------------------------------------------------------*/ \
  FOR_EACH_DIM(MK_DIM_UNITS, __VA_ARGS__) \
  \
  /*-----------------------------------------------------------------------*/ \
  /* Create a "DimLess" Type:                                              */ \
  /*-----------------------------------------------------------------------*/ \
  using DimLess = DimTypes::DimQ<0, 0, RepT, DimQ_MaxDims>; \
  \
  /*-----------------------------------------------------------------------*/ \
  /* Finally, generate a "DimQ" output function.                           */ \
  /* It uses the above-generated templates:                                */ \
  /*-----------------------------------------------------------------------*/ \
  MK_DIMQ_OUTPUT

//===========================================================================//
// Helper Methods for "DECLARE_DIMS":                                        //
//===========================================================================//
//---------------------------------------------------------------------------//
// "GET_DIM_NAME_COMMA", "GET_DIM_NAME":                                     //
//---------------------------------------------------------------------------//
// This is an "Action" for use with "FOR_EACH_DIM".
// Arg: (DimDcl=(DimName, FundUnitName [, OtherUnitDcl...])):
//
#ifdef  GET_DIM_NAME_COMMA
#undef  GET_DIM_NAME_COMMA
#endif
// NB: "DimDcl" is already enclosed in ()s, so "GET_DIM_NAME" is applied by
// juxtaposition. Also note the COMMA separator:
#define GET_DIM_NAME_COMMA(DimDcl) GET_DIM_NAME DimDcl ,

#ifdef  GET_DIM_NAME
#undef  GET_DIM_NAME
#endif
#define GET_DIM_NAME(    DimName, _FundUnitName, ...)  DimName

//---------------------------------------------------------------------------//
// "GET_FUND_UNIT", "MK_FUND_UNIT_STR":                                      //
//---------------------------------------------------------------------------//
// Similar to "GET_DIM_NAME" above:
//
#ifdef  GET_FUND_UNIT
#undef  GET_FUND_UNIT
#endif
#define GET_FUND_UNIT(   _DimName,  FundUnitName, ...)  FundUnitName

#ifdef  MK_FUND_UNIT_STR
#undef  MK_FUND_UNIT_STR
#endif
#define MK_FUND_UNIT_STR(_DimName,  FundUnitName, ...) #FundUnitName

//---------------------------------------------------------------------------//
// "GET_OTHER_UNITS":                                                        //
//---------------------------------------------------------------------------//
// Extracting OtherUnitNames from
// DimDcl =(DimName,  FundUnitName [, OtherUnitDcl...]),
// UnitDcl=(UnitName, UnitValInFundUnits) :
//
#ifndef GET_OTHER_UNITS_DCLS
#undef  GET_OTHER_UNITS_DCLS
#endif
#define GET_OTHER_UNITS_DCLS(_DimName, _FundUnitName, ...) __VA_ARGS__

#ifdef  GET_OTHER_UNITS
#undef  GET_OTHER_UNITS
#endif
#define GET_OTHER_UNITS(DimName, _FundUnitName, ...) \
  FOR_EACH_OTHER_UNIT(GET_UNIT_NAME_COMMA, DimName, __VA_ARGS__)

//---------------------------------------------------------------------------//
// "GET_UNIT_NAME*", "GET_UNIT_VAL":                                         //
//---------------------------------------------------------------------------//
#ifdef  GET_UNIT_NAME_COMMA
#undef  GET_UNIT_NAME_COMMA
#endif
// NB: "UnitDcl" is already enclosed in ()s, so "GET_UNIT_NAME" is applied by
// juxtaposition. Also note the COMMA separator:
#define GET_UNIT_NAME_COMMA(_DimName, UnitDcl) \
  GET_UNIT_NAME UnitDcl ,

#ifdef  GET_UNIT_NAME
#undef  GET_UNIT_NAME
#endif
#define GET_UNIT_NAME( UnitName, _UnitVal)  UnitName

#ifdef  GET_UNIT_VAL
#undef  GET_UNIT_VAL
#endif
#define GET_UNIT_VAL( _UnitName,  UnitVal)  UnitVal

//---------------------------------------------------------------------------//
// "MK_DIM_UNITS":                                                           //
//---------------------------------------------------------------------------//
// This is also an "Action" for use with "FOR_EACH_DIM".
// Arg: (DimDcl=(DimName,  FundUnitName [, OtherUnitDcl...])),
//       UnitDcl=(UnitName, UnitValInFundUnits)  :
// That is, again, "DimDcl" is already enclosed in ()s:
//
#ifdef  MK_DIM_UNITS
#undef  MK_DIM_UNITS
#endif
#define MK_DIM_UNITS(DimDcl) \
  /*-----------------------------------------------------------------------*/ \
  /* Check that this Dim is within the limits:                             */ \
  /*-----------------------------------------------------------------------*/ \
  static_assert(unsigned(DimsE::GET_DIM_NAME DimDcl) < DimQ_MaxDims, \
                "Too many Dimensions: " \
                STRINGIFY_NAME(GET_DIM_NAME DimDcl)); \
  /*-----------------------------------------------------------------------*/ \
  /* "{Dim}UnitsE" Enum:                                                   */ \
  /*-----------------------------------------------------------------------*/ \
  enum class MK_UNITS_ENUM_NAME DimDcl: unsigned  \
  { \
    GET_FUND_UNIT   DimDcl = 0, /* FundUnit has Code=0                    */  \
    GET_OTHER_UNITS DimDcl      /* NB: DimDcl is already enclosed in ()s  */  \
  }; \
  /*-----------------------------------------------------------------------*/ \
  /* "UnitName", "UnitScale"... specialisations for Fund and Other Units:  */ \
  /*-----------------------------------------------------------------------*/ \
  MK_UNIT_IMPL( GET_DIM_NAME DimDcl,   GET_FUND_UNIT DimDcl, 1.0)     \
  \
  FOR_EACH_OTHER_UNIT(MK_ANOTHER_UNIT, GET_DIM_NAME  DimDcl,          \
                      GET_OTHER_UNITS_DCLS   DimDcl) \
  \
  /*-----------------------------------------------------------------------*/ \
  /* Convenience Type for this Dim (with implicit FundUnit):               */ \
  /*-----------------------------------------------------------------------*/ \
  static_assert \
    (unsigned(MK_UNITS_ENUM_NAME DimDcl::GET_FUND_UNIT DimDcl) == 0); \
  using GET_DIM_NAME DimDcl = \
     DimTypes::DimQ   \
    <DimQ_Encs::DimExp(unsigned(DimsE::GET_DIM_NAME DimDcl)),    \
     DimQ_Encs::MkUnit(unsigned(DimsE::GET_DIM_NAME DimDcl), 0), \
     DimQ_RepT,       \
     DimQ_MaxDims>;   \
  /*-----------------------------------------------------------------------*/ \
  /* Simplified Conversion: "To_{DimName}" (assuming the FundUnit):        */ \
  /*-----------------------------------------------------------------------*/ \
  template<uint64_t E, uint64_t U> \
  constexpr DimTypes::DimQ   \
    <E, \
     DimQ_Encs::SetUnit(U, unsigned(DimsE:: GET_DIM_NAME DimDcl), 0),      \
     DimQ_RepT,        \
     DimQ_MaxDims>     \
  MK_CONV_FUNC_NAME0 DimDcl  \
  (DimTypes::DimQ<E, U, DimQ_RepT, DimQ_MaxDims> a_dimq) \
  { \
    constexpr unsigned Dim       = unsigned(DimsE::GET_DIM_NAME  DimDcl);  \
    constexpr unsigned OldUnit   = DimQ_Encs::GetFld(U, Dim);      \
    constexpr auto     ExpNumDen = \
      DimQ_Encs::GetNumerAndDenom (DimQ_Encs::GetFld(E, Dim));     \
    constexpr int      Numer     = ExpNumDen.first;  \
    constexpr unsigned Denom     = ExpNumDen.second; \
    return \
      DimTypes::DimQ \
      <E, DimQ_Encs::SetUnit(U, Dim, 0), DimQ_RepT, DimQ_MaxDims>  \
      ( \
        a_dimq.Magnitude() * \
        DimQ_Encs::FracPow<Numer, Denom> \
          /* OldScale / NewScale: */     \
          (UnitScale<Dim, OldUnit> / UnitScale<Dim, 0>) \
      ); \
  } \
  /*-----------------------------------------------------------------------*/ \
  /* For Convenience: Checking if "DimQ" is an "Elementary" Dim;           */ \
  /* The Unit Does Not Matter:                                             */ \
  /*-----------------------------------------------------------------------*/ \
  template<typename T> \
  constexpr inline bool MK_IS_ANY DimDcl = false;     \
  \
  template<uint64_t E, uint64_t U>                    \
  constexpr inline bool MK_IS_ANY DimDcl    \
    <DimTypes::DimQ<E, U, DimQ_RepT, DimQ_MaxDims>> = \
    (E == DimTypes::Bits::Encodings<DimQ_RepT, DimQ_MaxDims>::DimExp \
          (unsigned(DimsE::GET_DIM_NAME DimDcl)));

//---------------------------------------------------------------------------//
// "MK_ANOTHER_UNIT", "MK_UNIT_IMPL":                                        //
//---------------------------------------------------------------------------//
#ifdef  MK_ANOTHER_UNIT
#undef  MK_ANOTHER_UNIT
#endif
#define MK_ANOTHER_UNIT(DimName, UnitDcl)  \
        MK_UNIT_IMPL(   DimName, \
                        GET_UNIT_NAME  UnitDcl, GET_UNIT_VAL UnitDcl)

#ifdef  MK_UNIT_IMPL
#undef  MK_UNIT_IMPL
#endif
#define MK_UNIT_IMPL(DimName, UnitName, UnitVal) \
  static_assert(unsigned(MK_UNITS_ENTRY_NAME(DimName,UnitName)) <= \
                DimQ_Encs::PMask, \
                "Too many Units for " STRINGIFY_NAME(DimName)); \
  /*-----------------------------------------------------------------------*/ \
  /* "UnitNameStr" Specialisation for this Unit:                           */ \
  /*-----------------------------------------------------------------------*/ \
  template<> \
  inline constexpr char UnitNameStr   \
    <unsigned(DimsE::DimName), \
     unsigned(MK_UNITS_ENTRY_NAME(DimName,UnitName))>[] = \
    STRINGIFY_NAME(UnitName);  \
  /*-----------------------------------------------------------------------*/ \
  /* "UnitScale" Specialisation for this Unit:                             */ \
  /*-----------------------------------------------------------------------*/ \
  template<> \
  inline constexpr DimQ_RepT UnitScale     \
    <unsigned(DimsE::DimName), \
     unsigned(MK_UNITS_ENTRY_NAME(DimName,UnitName))> = DimQ_RepT(UnitVal);   \
  /*-----------------------------------------------------------------------*/ \
  /* Convenience Type and Literals for this Dim and Unit:                  */ \
  /*-----------------------------------------------------------------------*/ \
  using MK_DIMQ_TYPE_ALIAS(DimName, UnitName) = \
     DimTypes::DimQ \
    <DimQ_Encs::DimExp(unsigned(DimsE::DimName)), \
     DimQ_Encs::MkUnit(unsigned(DimsE::DimName),  \
                       unsigned(MK_UNITS_ENTRY_NAME(DimName,UnitName))), \
     DimQ_RepT,     \
     DimQ_MaxDims>; \
  /* NB: Literal suffix is the UnitName. Thus, identical UnitNames for     */ \
  /* different Units are NOT allowed (the stopper is precisely HERE!).     */ \
  /* This should normally be OK:                                           */ \
  constexpr MK_DIMQ_TYPE_ALIAS(DimName, UnitName) operator            \
    MK_LITERAL_OP(UnitName) (long double a_v)          \
    { return MK_DIMQ_TYPE_ALIAS(DimName, UnitName)(DimQ_RepT(a_v)); } \
  /*-----------------------------------------------------------------------*/ \
  /* Conversion function for DimQs: "To_{DimName}_{UnitName}:              */ \
  /*-----------------------------------------------------------------------*/ \
  template<uint64_t E, uint64_t U> \
  constexpr DimTypes::DimQ   \
    <E, \
     /* The Target Unit:  */ \
     DimQ_Encs::SetUnit        \
      (U, unsigned(DimsE::DimName), \
          unsigned(MK_UNITS_ENTRY_NAME(DimName,UnitName))),    \
     DimQ_RepT,     \
     DimQ_MaxDims>  \
  MK_CONV_FUNC_NAME(DimName,UnitName) \
  (DimTypes::DimQ<E, U, DimQ_RepT, DimQ_MaxDims> a_dimq)       \
  { \
    constexpr unsigned Dim       = unsigned(DimsE::DimName);   \
    constexpr unsigned OldUnit   = DimQ_Encs::GetFld(U, Dim);  \
    constexpr unsigned NewUnit   = \
      unsigned(MK_UNITS_ENTRY_NAME(DimName,UnitName));         \
    constexpr auto     ExpNumDen = \
      DimQ_Encs::GetNumerAndDenom(DimQ_Encs::GetFld(E, Dim));  \
    constexpr int      Numer     = ExpNumDen.first;            \
    constexpr unsigned Denom     = ExpNumDen.second;           \
    return \
      DimTypes::DimQ \
      <E, DimQ_Encs::SetUnit(U, Dim, NewUnit), DimQ_RepT, DimQ_MaxDims> \
      ( \
        a_dimq.Magnitude() * \
        DimQ_Encs::FracPow<Numer, Denom>                       \
          /* OldScale / NewScale: */                           \
          (UnitScale<Dim, OldUnit> / UnitScale<Dim, NewUnit>)  \
      ); \
  }

//===========================================================================//
// NAME CONCATENATIONS AND STRINGIFICATION:                                  //
//===========================================================================//
// XXX: Why are they required in the first place?   That is, why can't we use
// "##" concatenations directly in the calling macros? -- For some reason, in
// the latter case,  component names may not get fully evaluated,  so we need
// yet another level of macro invocation to force evaluation:
//
#ifdef  MK_UNITS_ENUM_NAME
#undef  MK_UNITS_ENUM_NAME
#endif
#define MK_UNITS_ENUM_NAME( DimName, _FundUnitName, ...) DimName##UnitsE

#ifdef  MK_UNITS_ENTRY_NAME
#undef  MK_UNITS_ENTRY_NAME
#endif
#define MK_UNITS_ENTRY_NAME(DimName, UnitName)  DimName##UnitsE::UnitName

#ifdef  MK_DIMQ_TYPE_ALIAS
#undef  MK_DIMQ_TYPE_ALIAS
#endif
#define MK_DIMQ_TYPE_ALIAS(DimName, UnitName)   DimName##_##UnitName

#ifdef  MK_IS_ANY
#undef  MK_IS_ANY
#endif
#define MK_IS_ANY(DimName, _FundUnitName, ...)  IsAny##DimName

#ifdef  MK_LITERAL_OP
#undef  MK_LITERAL_OP
#endif
#define MK_LITERAL_OP(UnitName)                 ""_##UnitName

#ifdef  MK_CONV_FUNC_NAME
#undef  MK_CONV_FUNC_NAME
#endif
#define MK_CONV_FUNC_NAME(DimName, UnitName)    To_##DimName##_##UnitName

#ifdef  MK_CONV_FUNC_NAME0
#undef  MK_CONV_FUNC_NAME0
#endif
#define MK_CONV_FUNC_NAME0(DimName, ...)        To_##DimName

#ifdef  STRINGIFY_NAME
#undef  STRINGIFY_NAME
#endif
#define STRINGIFY_NAME( SomeName) STRINGIFY_NAME1(SomeName)

#ifdef  STRINGFFY_NAME1
#undef  STRINGIFY_NAME1
#endif
#define STRINGIFY_NAME1(SomeName) #SomeName

//===========================================================================//
// "DimQ" Output Function. NB: it is NOT a "constexpr"!                      //
//===========================================================================//
//---------------------------------------------------------------------------//
// "MK_DIMQ_OUTPUT":                                                         //
//---------------------------------------------------------------------------//
#ifdef  MK_DIMQ_OUTPUT
#undef  MK_DIMQ_OUTPUT
#endif
#define MK_DIMQ_OUTPUT \
  /*-----------------------------------------------------------------------*/ \
  /* "Put": Safe outputting of a "DimQ" into a given Buffer:               */ \
  /*-----------------------------------------------------------------------*/ \
  template<uint64_t E, uint64_t U>  \
  char* Put \
  ( \
    DimTypes::DimQ<E, U, DimQ_RepT, DimQ_MaxDims> a_dimq, \
    char*                                         a_buff, \
    char const*                                   a_end   \
  ) \
  { \
    char* curr = a_buff;           \
    /* Available buffer size: */   \
    int   N    = 0;                \
    CHECK_N /* Sets N */           \
    /* First, the Magnitude (Real or Complex): */ \
    curr = DimQ_Encs::PutMagnitude(curr, N, a_dimq.Magnitude());         \
    CHECK_N \
    /* Now the Units and Exponents, up to the AbsoluteMaxDims-1 = 8: */  \
    /* Dims >= DimQ_MaxDims will be safely ignored:                  */  \
    PUT_UNIT_AND_EXP_STR(0) \
    PUT_UNIT_AND_EXP_STR(1) \
    PUT_UNIT_AND_EXP_STR(2) \
    PUT_UNIT_AND_EXP_STR(3) \
    PUT_UNIT_AND_EXP_STR(4) \
    PUT_UNIT_AND_EXP_STR(5) \
    PUT_UNIT_AND_EXP_STR(6) \
    PUT_UNIT_AND_EXP_STR(7) \
    PUT_UNIT_AND_EXP_STR(8) \
    /* For safety, always 0-terminate the buffer: */ \
    CHECK_N \
    *curr = '\0'; \
    return  curr; \
  } \
  /*-----------------------------------------------------------------------*/ \
  /* "ToStr": Conversion of "DimQ" into a fixed-size string:               */ \
  /*-----------------------------------------------------------------------*/ \
  /* XXX: Are 64 bytes OK for most use cases? */     \
  template<size_t Sz = 64, uint64_t E, uint64_t U>   \
  inline std::array<char, Sz> ToStr   \
    (DimTypes::DimQ<E, U, DimQ_RepT, DimQ_MaxDims> a_dimq) \
  { \
    std::array<char, Sz>   buff;        \
    char const*            buffEnd = buff.data() + Sz; \
    [[maybe_unused]] char* realEnd = Put<E, U>(a_dimq, buff.data(), buffEnd); \
    assert(realEnd < buffEnd && *realEnd == '\0');   \
    return buff;    \
  } \
  /*-----------------------------------------------------------------------*/ \
  /* Output: "operator<<":                                                 */ \
  /*-----------------------------------------------------------------------*/ \
  template<uint64_t E, uint64_t U> \
  std::ostream& operator<< \
  ( \
     std::ostream& a_os,   \
     DimTypes::DimQ<E, U, DimQ_RepT, DimQ_MaxDims> a_dimq \
  ) \
  { return a_os << ToStr(a_dimq).data(); }

//---------------------------------------------------------------------------//
// "PUT_UNIT_AND_EXP_STR":                                                   //
//---------------------------------------------------------------------------//
// For use in the context of "MK_DIMQ_OUTPUT":
//
#ifdef  PUT_UNIT_AND_EXP_STR
#undef  PUT_UNIT_AND_EXP_STR
#endif
#define PUT_UNIT_AND_EXP_STR(Dim) \
  /* NB: We may have Dim >= DimQ_MaxDims, but in that case, we do nothing: */ \
  if constexpr(Dim < DimQ_MaxDims) \
  { \
    constexpr auto     NumDen = DimQ_Encs::GetNumerAndDenom \
                               (DimQ_Encs::GetFld(E,Dim));  \
    constexpr int      Numer  = NumDen.first;   \
    constexpr unsigned Denom  = NumDen.second;  \
    /* NB: Units are only relevant if the exponent is non-0: */ \
    if constexpr(Numer != 0) \
    { \
      /* Output the UnitName, preceded by a space: */ \
      curr += snprintf(curr, size_t(N), " %s",   \
                       UnitNameStr<Dim, DimQ_Encs::GetFld(U,Dim)>); \
      CHECK_N \
      /* Output the Exponent, unless it is 1: */ \
      if constexpr(!(Numer == 1 && Denom == 1))  \
      { \
        constexpr bool brackets =  Numer < 0 || Denom != 1;  \
        curr += snprintf(curr, size_t(N), brackets ? "^(%d" : "^%d", Numer); \
        CHECK_N \
        if constexpr(Denom != 1) \
        { \
          curr += snprintf(curr, size_t(N), "/%u", Denom);   \
          CHECK_N \
        } \
        if constexpr(brackets)   \
        { \
          *curr = ')'; \
          ++curr; \
          CHECK_N \
        } \
      } \
    } \
  }
//---------------------------------------------------------------------------//
// "CHECK_N":                                                                //
//---------------------------------------------------------------------------//
#ifdef  CHECK_N
#undef  CHECK_N
#endif
#define CHECK_N \
  N = int(a_end - curr);     \
  /* Always allow a reserve of at least 1 byte: */ \
  if (N <= 1) \
    throw std::runtime_error("Put(DimQ): Buffer OverFlow");

//===========================================================================//
// "MK_DIMS_FMT":                                                            //
//===========================================================================//
// Creates a Formatter (for use with "std::format") for "DimQ"s previously dec-
// lared by "DECLARE_DIMS". Because the Formatter requires specialisation of a
// class declared in "std", it cannot in general be created  in "DECLARE_DIMS"
// itself, since the latter macro may be invoked inside any namespace "DIMS_NS".
// We thus need a separate macro, "MK_DIMS_FMT". It is to be invoked OUTSIDE of
// any namespaces, and after the  "DECLARE_DIMS" declaration.   The name of the
// namespace ("DIMS_NS") in which "DECLARE_DIMS" was invoked, if exists, must be
// passed to "MK_DIMS_FMT". This is an unpleasant hack...
//
#ifdef  MK_DIMS_FMT
#undef  MK_DIMS_FMT
#endif
#define MK_DIMS_FMT(DIMS_NS) \
  namespace std \
  { \
    template<uint64_t E, uint64_t U> \
    struct formatter<DimTypes::DimQ  \
      <E, U, DIMS_NS::DimQ_RepT, DIMS_NS::DimQ_MaxDims>>     \
    { \
      constexpr auto parse(std::format_parse_context& a_ctx) \
        { return a_ctx.begin(); } \
      \
      auto format(DimTypes::DimQ  \
                  <E, U, DIMS_NS::DimQ_RepT, DIMS_NS::DimQ_MaxDims> a_dimq, \
                  format_context&                                   a_ctx)  \
      const \
      { \
        return std::format_to(a_ctx.out(), "{}", \
               DIMS_NS::ToStr(a_dimq).data());   \
      } \
    };  \
  }

//===========================================================================//
// Run-Time Checks:                                                          //
//===========================================================================//
//---------------------------------------------------------------------------//
// "CHECK_ONLY":                                                             //
//---------------------------------------------------------------------------//
// "CHECK_ONLY": Symbol is used unless UNCHECKED_MODE (and thus NDEBUG) are set:
// To make sure: UNCHECKED_MODE => NDEBUG, but not other way round:
//
#ifdef CHECK_ONLY
#undef CHECK_ONLY
#endif
#if UNCHECKED_MODE
#  define CHECK_ONLY(x)           // UnChecked Mode only
#  ifndef NDEBUG
#  define NDEBUG 1
#  endif
#else
#  define CHECK_ONLY(x) x         // Release, RelWithDebInfo, Debug
#endif

//---------------------------------------------------------------------------//
// "DEBUG_ONLY":                                                             //
//---------------------------------------------------------------------------//
// "DEBUG_ONLY": Symbol is used unless NDEBUG is set, which is a STRONGER cond
// compared to "CHECK_ONLY":
#ifdef DEBUG_ONLY
#undef DEBUG_ONLY
#endif
#ifdef NDEBUG
#  define DEBUG_ONLY(x)           // UnChecked or Release Mode
#else
#  define DEBUG_ONLY(x) x         // RelWithDebInfo or Debug
#endif

//---------------------------------------------------------------------------//
// "UNUSED_PARAM" (unconditionally):                                         //
//---------------------------------------------------------------------------//
#ifdef  UNUSED_PARAM
#undef  UNUSED_PARAM
#endif
#define UNUSED_PARAM(x)

//---------------------------------------------------------------------------//
// Branch Prediction:                                                        //
//---------------------------------------------------------------------------//
// NB: "LIKELY" and "UNLINKEY" macros are to distinguish deeply-generic and
// deeply-degenerate cases, so the probabilities of each case are set high.
// Indended to be used at run-time only,  so do not use them in "consnexpr"
// functions:
//
#ifndef  LIKELY
# define LIKELY(expr)   __builtin_expect_with_probability((expr), 1, 0.999)
#endif
#ifndef  UNLIKELY
# define UNLIKELY(expr) __builtin_expect_with_probability((expr), 0, 0.999)
#endif

