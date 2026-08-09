// vim:ts=2:et
//===========================================================================//
//                          "DimTypes/DimTypes.hpp":                         //
//          Types for Physical Dimensions and Dimensioned Quantities         //
//===========================================================================//
#pragma  once
#include "ValidRepTs.h"
#include "CEMaths.hpp"
#include "Encodings.hpp"
#include "CompilerMacros.h"
#include "DeclareDims.h"

namespace DimTypes
{
  // OUR STANDARD LOOKUP / OVERLOADING PATTERN (similar to "CEMaths.hpp"):
  // Bring in an entity from CEMaths, but overloaded (same-named) entities may
  // also exist in the context -- overloading resolution will then be performed:
  //
  using CEMaths::IsComplex;
  using CEMaths::DefaultTol;
  using CEMaths::ApproxEqual;

  //=========================================================================//
  // "DimQ": The External Dimensioned Type:                                  //
  //=========================================================================//
  // "DimQ" is actually a Field "T" statically parameterised by the encoded
  // Dimension Exponent "E" and the vector of units "U":
  // NB: "RepT" and "DimQ" are everywhere passed BY COPY (except for result vals
  // of in-place update operators),  which is probably optimal even if "RepT" is
  // a complex or multiple-precision type, to provide better data locality;
  // "MaxDims" can be 7, 8 or 9:
  //
  template<uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  requires IsValidDimQRepT<RepT>
  class DimQ
  {
  private:
    RepT  m_val; // The actual value (magnitude)

    // Encodings to be used:
    using En =   Bits::Encodings<RepT, MaxDims>;

  public:
    //=======================================================================//
    // Ctors and Accessors:                                                  //
    //=======================================================================//
    // Default Ctor:  Initialisation to 0:
    constexpr DimQ(): m_val(0.0) {}

    // Non-Default Ctor:
    // Lifting a value of type "R" into "DimQ".  We DELIBERATELY DO NOT INFER
    // the "RepT" from the type of the Ctor arg, as this may cause ambiguities
    // in some cases.  "R" is required to be a "valid RepT" (but not nrcessar-
    // ili the "RepT" itself) to avoid a conflict with the Copy Ctor below
    // (where "R" is "DimQ" itself):
    template<typename  R> requires(IsValidDimQRepT<R>)
    constexpr explicit DimQ(R const& val)
    : m_val(RepT(val)) {}

    // Copy Ctor and assignment are auto-generated:
    constexpr DimQ            (DimQ const&) = default;
    constexpr DimQ& operator= (DimQ const&) = default;

    // If the assignment LHS is DimLess, the RHS may be of any floating-point
    // type "R"   (again, this constraint is required to avoid confusion with
    // the "DimQ" itself):
    template<typename R> requires(IsValidDimQRepT<R> && E == 0)
    constexpr DimQ& operator= (R const& a_right)
    {
      m_val = RepT(a_right);
      return *this;
    }

    // Non-Default Ctor:
    // Conversion FROM another Rep (but E, U and MaxDims must match); this is
    // actually an "Extended Copy Ctor":
    template<typename  FromR> requires(IsValidDimQRepT<FromR>)
    constexpr explicit DimQ(DimQ<E, U, FromR, MaxDims> const& a_right)
      : m_val(RepT(a_right.Magnitude()))
      {}

    // The dimensioned unit for the given dimensioned quantity. NB: This method
    // is not a ctor: it takes an existing "DimQ" and returns a similar one but
    // with the unitary magnitude:
    //
    constexpr DimQ     UnitOf()       const { return DimQ(1.0); }

    // Getting the exponent code (this is primarily for testing):
    constexpr uint64_t GetDimsCode () const { return E; }

    // Getting the Units code (also for testing):
    constexpr uint64_t GetUnitsCode() const { return U; }

    // The Magnitude of the dimensioned quantity (ie its value expressed in the
    // corresp dimensioned units). The magnitude is by itself dimension-less.
    // XXX: Direct access to the dimension-less magnitude  may seem  to be an
    // unsafe feature but it could be emulated anyway (by dividing the qty by
    // its unit), so we provide the direct access for efficiency:
    //
    constexpr RepT const& Magnitude() const { return m_val; }

    // The following function converts "DimQ" to the Underlying Field Type
    // ("RepT"), but ONLY if the arg is actually dimension-less:
    template<typename R> requires(IsValidDimQRepT<R> && E == 0)
    constexpr explicit   operator R() const
      { return R(m_val); }
    // NB:  In addition, for some types (eg Angles expressed in Radians), the
    // user may allow direct conversion into "RepT" by specialising the above
    // operator method...

    // Creating an Invalid value:
    constexpr static DimQ Invalid() { return DimQ(NAN); }

    //=======================================================================//
    // Arithmetic Operations -- Dimensions Unchanged:                        //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Addition and Subtraction of "DimQ"s:                                  //
    //-----------------------------------------------------------------------//
    // This is only possible if the dimensions and RepT are same ("E"), and the
    // units ("U" and "V") can be unified  (ie same units for non-0 exponents),
    // in which case the 1st arg's units can still be used:
    //-----------------------------------------------------------------------//
    // Addition:                                                             //
    //-----------------------------------------------------------------------//
    template<uint64_t V>
    requires(En::UnitsOK(E,U,V))
    constexpr  DimQ operator+ (DimQ<E,V,RepT,MaxDims> const& a_right) const
      { return DimQ(m_val   + a_right.Magnitude()); }

    template<uint64_t V>
    requires(En::UnitsOK(E,U,V))
    constexpr DimQ& operator+= (DimQ<E,V,RepT,MaxDims> const& a_right)
    {
      m_val += a_right.Magnitude();
      return *this;
    }

    // Addition of a "scalar" is permitted only for a dimension-less qty:
    template<typename R>
    requires(IsValidDimQRepT<R> && E == 0)
    constexpr DimQ operator+(R const& a_right) const
    { return  DimQ(m_val   + RepT(a_right)); }

    template<typename R>
    requires(IsValidDimQRepT<R> && E == 0)
    constexpr DimQ& operator+=(R const& a_right)
    {
      m_val += RepT(a_right);
      return *this;
    }

    //-----------------------------------------------------------------------//
    // Subtraction of "DimQs":                                               //
    //-----------------------------------------------------------------------//
    // Same constraints as for Addition:
    //
    template<uint64_t V>
    requires(En::UnitsOK(E,U,V))
    constexpr  DimQ  operator- (DimQ<E,V,RepT,MaxDims> const& a_right) const
      { return DimQ(m_val - a_right.Magnitude()); }

    template<uint64_t V>
    requires(En::UnitsOK(E,U,V))
    constexpr  DimQ& operator-=(DimQ<E,V,RepT,MaxDims> const& a_right)
    {
      m_val -= a_right.Magnitude();
      return *this;
    }

    template<typename R>
    requires(IsValidDimQRepT<R> && E == 0)
    constexpr DimQ operator-(R const& a_right) const
    { return  DimQ(m_val - RepT(a_right)); }

    template<typename R>
    requires(IsValidDimQRepT<R> && E == 0)
    constexpr DimQ& operator-=(R const& a_right)
    {
      m_val -= RepT(a_right);
      return *this;
    }

    //-----------------------------------------------------------------------//
    // Multiplication and Division by some "R":                              //
    //-----------------------------------------------------------------------//
    // NB: the following method can be used to create new dimensioned values
    // from the fundamental units. XXX: when used this way, it is slighly in-
    // efficient because of multiplication by 1 and extra copying;   but the
    // impact of that should be negligible:
    //
    // Multiplication:
    template<typename R> requires(IsValidDimQRepT<R>)
    constexpr  DimQ operator*    (R const& a_right) const
      { return DimQ(m_val   * RepT(a_right)); }

    template<typename R> requires(IsValidDimQRepT<R>)
    constexpr  DimQ& operator*=  (R const& a_right)
    {
      m_val *= RepT(a_right);
      return  *this;
    }

    // In particular, Unary Negation:
    constexpr  DimQ operator- () const
      { return DimQ(- m_val); }

    // Division by "R":
    template<typename R> requires(IsValidDimQRepT<R>)
    constexpr  DimQ operator/    (R const& a_right) const
      { return DimQ(m_val   / RepT(a_right)); }

    template<typename R> requires(IsValidDimQRepT<R>)
    constexpr  DimQ& operator/=  (R const& a_right)
    {
      m_val /= RepT(a_right);
      return *this;
    }

    //=======================================================================//
    // Arithmetic operations which result in Dimensions change:              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Multiplication of "DimQ"s:                                            //
    //-----------------------------------------------------------------------//
    // Dimension exponents are added (or subtracted) up; units must unify;
    // afterwards, units of zero-exp dims are reset, otherwise false type-check-
    // ing failures may occur:
    //
    template<uint64_t F, uint64_t V>
    constexpr DimQ<En::AddExp(E,F),
                   En::CleanUpUnits(En::AddExp(E,F), En::UnifyUnits(E,F,U,V)),
                   RepT,
                   MaxDims>
    operator* (DimQ<F, V, RepT, MaxDims> const& a_right) const
    {
      return  DimQ<En::AddExp(E,F),
                   En::CleanUpUnits(En::AddExp(E,F), En::UnifyUnits(E,F,U,V)),
                   RepT,
                   MaxDims>
             (m_val * a_right.Magnitude());
    }

    //-----------------------------------------------------------------------//
    // Division of "DimQ"s:                                                  //
    //-----------------------------------------------------------------------//
    // Dimension exponents are subtracted. Same treatment of units as for mult:
    template<uint64_t F, uint64_t V>
    constexpr DimQ<En::SubExp(E,F),
                   En::CleanUpUnits(En::SubExp(E,F), En::UnifyUnits(E,F,U,V)),
                   RepT,
                   MaxDims>
    operator/ (DimQ<F, V, RepT, MaxDims> const& a_right) const
    {
      return  DimQ<En::SubExp(E,F),
                   En::CleanUpUnits(En::SubExp(E,F), En::UnifyUnits(E,F,U,V)),
                   RepT,
                   MaxDims>
             (m_val / a_right.Magnitude());
    }

    //-----------------------------------------------------------------------//
    // Comparison operators:                                                 //
    //-----------------------------------------------------------------------//
    // Require same representation, same dimensions and unifyable units.
    // NB: different requirements for DIMQ_EQ_CMP and DIMQ_ORD_CMP:
    //
#   ifdef  DIMQ_EQ_CMP
#   undef  DIMQ_EQ_CMP
#   endif
#   define DIMQ_EQ_CMP(Op) \
    template<uint64_t V>   \
    requires(En::UnitsOK(E,U,V)) \
    constexpr bool operator Op(DimQ<E, V, RepT, MaxDims> const& a_right) const \
    { return m_val Op a_right.Magnitude(); }

    DIMQ_EQ_CMP(==)
    DIMQ_EQ_CMP(!=)
#   undef  DIMQ_EQ_CMP

#   ifdef  DIMQ_ORD_CMP
#   undef  DIMQ_ORD_CMP
#   endif
#   define DIMQ_ORD_CMP(Op) \
    template<uint64_t V>   \
    requires(En::UnitsOK(E,U,V) && !IsComplex<RepT>) \
    constexpr bool operator Op(DimQ<E, V, RepT, MaxDims> const& a_right) const \
    { return m_val Op a_right.Magnitude(); }

    DIMQ_ORD_CMP(>)
    DIMQ_ORD_CMP(>=)
    DIMQ_ORD_CMP(<)
    DIMQ_ORD_CMP(<=)
#   undef  DIMQ_ORD_CMP

    //-----------------------------------------------------------------------//
    // Approximate Equality:                                                 //
    //-----------------------------------------------------------------------//
    template<uint64_t V>
    requires(En::UnitsOK(E,U,V) && !IsComplex<RepT>)
    constexpr bool ApproxEquals
    (
      DimQ<E, V, RepT, MaxDims> const& a_right,
      RepT                      const& a_tol = DefaultTol<RepT>
    )
    const
    { return ApproxEqual(m_val, a_right.Magnitude(), a_tol); }
  };
  // End "DimQ" class

  //=========================================================================//
  // NameSpace-Level Functions:                                              //
  //=========================================================================//
  // IMPORTANT: DO NOT make them "friends", as this would limit the template
  // substitution / overloading resolution ("friend hiding"):!
  //-------------------------------------------------------------------------//
  // Conversion TO another Rep (E, U and MaxDims are unchanged):             //
  //-------------------------------------------------------------------------//
  template<typename ToR,   // "ToR" must be explicit, other params are inferred
           uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr
     DimQ<E, U, ToR,  MaxDims> ToNewRepT
    (DimQ<E, U, RepT, MaxDims> const& a_right)
    // This is just a shortcut to the "Extended Copy Ctor":
    { return DimQ<E, U, ToR, MaxDims>(a_right); }

  //-------------------------------------------------------------------------//
  // "Abs", "UniformNorm", "Floor", "Ceil", "Round": "DimQ" -> "DimQ":       //
  //-------------------------------------------------------------------------//
# ifdef  DIMQ_ROUND_FUNC
# undef  DIMQ_ROUND_FUNC
# endif
# define DIMQ_ROUND_FUNC(FuncName)   \
  /* Bring in CMaths::FuncName (though FuncName may be overloaded): */  \
  using CEMaths::FuncName;            \
  \
  template<uint64_t E,  uint64_t  U, typename RepT, unsigned MaxDims>   \
  constexpr DimQ<E, U, RepT, MaxDims> \
  FuncName( DimQ<E, U, RepT, MaxDims>  const& a_right) \
  requires(!IsComplex<RepT>) \
    { return DimQ<E, U, RepT, MaxDims>(FuncName(a_right.Magnitude()));  }

  DIMQ_ROUND_FUNC(Abs)
  DIMQ_ROUND_FUNC(UniformNorm)
  DIMQ_ROUND_FUNC(Floor)
  DIMQ_ROUND_FUNC(Ceil)
  DIMQ_ROUND_FUNC(Round)
# undef  DIMQ_ROUND_FUNC

  //-------------------------------------------------------------------------//
  // "ApproxEqual" (w/o "s", as in "CEMaths"):                               //
  //-------------------------------------------------------------------------//
  template<uint64_t E, uint64_t U, uint64_t V, typename RepT, unsigned MaxDims>
  requires(Bits::Encodings<RepT, MaxDims>::UnitsOK(E, U, V) && !IsComplex<RepT>)
  constexpr bool ApproxEqual
  (
    DimQ<E, U, RepT, MaxDims> const& a_left,
    DimQ<E, V, RepT, MaxDims> const& a_right,
    RepT                      const& a_tol = DefaultTol<RepT>
  )
  { return a_left.ApproxEquals(a_right, a_tol); }

  //-------------------------------------------------------------------------//
  // Addition, Subraction:                                                   //
  //-------------------------------------------------------------------------//
  // XXX: For "a_right", E==0 is required; U != 0 is allowed:
  //
  template<uint64_t U, typename R, typename RepT, unsigned MaxDims>
  requires(IsValidDimQRepT<R>)
  constexpr DimQ<0, U, RepT, MaxDims>
  operator+
  (
    R const&                           a_left,
    DimQ<0, U, RepT, MaxDims> const&   a_right
  )
  { return a_right + a_left; }

  template<uint64_t U, typename R, typename RepT, unsigned MaxDims>
  requires(IsValidDimQRepT<R>)
  constexpr DimQ<0, U, RepT, MaxDims>
  operator-
  (
    R const&                           a_left,
    DimQ<0, U, RepT, MaxDims> const&   a_right
  )
  {
    return DimQ<0, U, RepT, MaxDims>
           (RepT(a_left) - a_right.Magnitude());
  }

  //-------------------------------------------------------------------------//
  // Multiplication, Division:                                               //
  //-------------------------------------------------------------------------//
  template<typename R, uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  requires(IsValidDimQRepT<R>)
  constexpr
  DimQ<E, U, RepT, MaxDims> operator*
    (R const& a_left, DimQ<E, U, RepT, MaxDims> const& a_right)
    { return DimQ<E, U, RepT, MaxDims>(RepT(a_left)  * a_right.Magnitude()); }

  template<typename R, uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  requires(IsValidDimQRepT<R>)
  constexpr auto    // Using "auto" here for brievity
  operator/ (R const& a_left, DimQ<E, U,  RepT, MaxDims> const& a_right)
  {
    using  En = Bits::Encodings<RepT,  MaxDims>;
    return DimQ<En::SubExp(0UL,E),  U, RepT, MaxDims>
               (RepT(a_left) / a_right.Magnitude());
  }

  //-------------------------------------------------------------------------//
  // Integral and Rational Powers:                                           //
  //-------------------------------------------------------------------------//
  // The power must be given by constant exprs, so it is analysable at compile
  // time:
  // NB: "CleanUpUnits" is required in case of M==0:
  //-------------------------------------------------------------------------//
  // "IPow": Integral Power:                                                 //
  //-------------------------------------------------------------------------//
  template<int M,   // "M" must be explicit, other params are inferred:
           uint64_t E, uint64_t  U, typename RepT, unsigned MaxDims>
  constexpr auto
  IPow(DimQ<E, U, RepT, MaxDims> const& a_right)
  {
    using  En = Bits::Encodings<RepT,  MaxDims>;
    return
      DimQ<En::MultExp(E,M), En::CleanUpUnits(En::MultExp(E,M), U),
           RepT,             MaxDims>
          (En::template IntPow<M>(a_right.Magnitude()));
  }

  //-------------------------------------------------------------------------//
  // "Sqr" and "Cube" are particularly-important cases of "IPow":            //
  //-------------------------------------------------------------------------//
  template<uint64_t  E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto
  Sqr(DimQ<E, U, RepT, MaxDims> const& a_right)
    { return IPow<2>(a_right); }

  template<uint64_t  E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto
  Cube(DimQ<E, U, RepT, MaxDims> const& a_right)
    { return IPow<3>(a_right); }

  // Bring in CEMaths::{Sqr,Cube} (though they may be overloaded):
  using CEMaths::Sqr;
  using CEMaths::Cube;

  //-------------------------------------------------------------------------//
  // "RPow": (General) Rational Power:                                       //
  //-------------------------------------------------------------------------//
  // XXX: In general it is NOT "constexpr" for C++ < 26, but it is "constexpr"
  // if it can be reduced to a composition of "SqRt" and "CbRt":
  //
  template<int M, int N, // "M" and "N" must be explicit, others are inferred
           uint64_t   E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto
  RPow(DimQ<E, U, RepT, MaxDims> const& a_right)
  {
    // Multiples of "PMod" (including 0) are uninvertible in the Z_{PMod} rep:
    using  En = Bits::Encodings<RepT,  MaxDims>;
    static_assert(N > 0 && N % En::IPMod != 0, "RPow: Invalid Denom");
    return
      DimQ<En::DivExp(En::MultExp(E,M),N),
           En::CleanUpUnits(En::DivExp(En::MultExp(E,M),N), U),
           RepT,
           MaxDims>
          (En::template FracPow<M, N>(a_right.Magnitude()));
  }

  //-------------------------------------------------------------------------//
  // "SqRt" and "CbRt": "RPow" Shortcuts:                                    //
  //-------------------------------------------------------------------------//
  // Here M==1, so "CleanUpUnits" is not required (it is only needed for M==0):
  //
  template<uint64_t  E,  uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto
  SqRt(DimQ<E, U, RepT, MaxDims> const& a_right)
    { return RPow<1,2>(a_right); }

  template<uint64_t  E,  uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto
  CbRt(DimQ<E, U, RepT, MaxDims> const& a_right)
    { return RPow<1,3>(a_right); }

  // Bring in CEMaths::{SqRt,CbRt} (though they may be overloaded):
  using CEMaths::SqRt;
  using CEMaths::CbRt;

  //-------------------------------------------------------------------------//
  // Unary Mathematical Funcs: DimLess->RepT:                                //
  //-------------------------------------------------------------------------//
# ifdef  DIMLESS_UNARY_FUNC
# undef  DIMLESS_UNARY_FUNC
# endif
# define DIMLESS_UNARY_FUNC(FuncName)    \
  /* Bring in CEMaths::FuncName (though it may be overloaded): */      \
  using CEMaths::FuncName; \
  \
  /* We must have E == 0, but U != 0 is allowed as above:      */      \
  template<uint64_t U, typename RepT, unsigned MaxDims>    \
  constexpr RepT FuncName(DimQ<0, U, RepT, MaxDims> const& a_right)    \
  { \
    /* NB: The following call to "FuncName" MAY OR MAY NOT BE from */  \
    /* "CEMaths"! */ \
    return FuncName(a_right.Magnitude()); \
  }
  DIMLESS_UNARY_FUNC(Exp)
  DIMLESS_UNARY_FUNC(Log)
  DIMLESS_UNARY_FUNC(Log2)
  DIMLESS_UNARY_FUNC(Log10)
  DIMLESS_UNARY_FUNC(Cos)
  DIMLESS_UNARY_FUNC(Sin)
  DIMLESS_UNARY_FUNC(Tan)
  DIMLESS_UNARY_FUNC(ATan)
  DIMLESS_UNARY_FUNC(ASin)
  DIMLESS_UNARY_FUNC(ACos)
  DIMLESS_UNARY_FUNC(CosH)
  DIMLESS_UNARY_FUNC(SinH)
  DIMLESS_UNARY_FUNC(TanH)
  DIMLESS_UNARY_FUNC(ACosH)
  DIMLESS_UNARY_FUNC(ASinH)
  DIMLESS_UNARY_FUNC(ATanH)
# undef DIMLESS_UNARY_FUNC

  //-------------------------------------------------------------------------//
  // Predicates (XXX: currently for Real "RepT" only):                       //
  //-------------------------------------------------------------------------//
# ifdef  DIMQ_PRED
# undef  DIMQ_PRED
# endif
# define DIMQ_PRED(Pred) \
  /* Bring in CEMaths::Pred: */ \
  using CEMaths::Pred; \
  \
  template<uint64_t E, uint64_t U, typename RepT, unsigned MaxDims> \
  constexpr bool Pred(DimQ<E, U, RepT, MaxDims>   const& a_right)   \
  requires(!IsComplex<RepT>) \
  { return Pred(a_right.Magnitude()); }

  DIMQ_PRED(IsZero)
  DIMQ_PRED(IsFinite)
  DIMQ_PRED(IsInf)
  DIMQ_PRED(IsNaN)
  DIMQ_PRED(IsNeg)
  DIMQ_PRED(IsPos)
# undef  DIMQ_PRED

  //-------------------------------------------------------------------------//
  // "ATan2" on "DimQs":                                                     //
  //-------------------------------------------------------------------------//
  using CEMaths::ATan2;

  template<uint64_t E, uint64_t U, uint64_t V, typename RepT,
           unsigned MaxDims>
  requires(Bits::Encodings<RepT,MaxDims>::UnitsOK(E,U,V))
  constexpr  RepT   ATan2
    (DimQ<E, U, RepT, MaxDims> const& a_y,
     DimQ<E, V, RepT, MaxDims> const& a_x)
  requires(!IsComplex<RepT>)
  { return ATan2(a_y.Magnitude(), a_x.Magnitude()); }
}
// End namespace DimTypes
//===========================================================================//
// Deferred Import of "DimTypes" Elementary Functions (on "DimQ"):           //
//===========================================================================//
#ifdef  IMPORT_DIMTYPES_FUNCS
#undef  IMPORT_DIMTYPES_FUNCS
#endif
#define IMPORT_DIMTYPES_FUNCS           \
  using DimTypes::ToNewRepT;            \
  using DimTypes::Abs;                  \
  using DimTypes::Floor;                \
  using DimTypes::Ceil;                 \
  using DimTypes::Round;                \
  using DimTypes::ApproxEqual;          \
  using DimTypes::IPow;                 \
  using DimTypes::Sqr;                  \
  using DimTypes::Cube;                 \
  using DimTypes::RPow;                 \
  using DimTypes::SqRt;                 \
  using DimTypes::CbRt;                 \
  using DimTypes::Exp;                  \
  using DimTypes::Log;                  \
  using DimTypes::Log2;                 \
  using DimTypes::Log10;                \
  using DimTypes::Cos;                  \
  using DimTypes::Sin;                  \
  using DimTypes::Tan;                  \
  using DimTypes::ATan;                 \
  using DimTypes::ASin;                 \
  using DimTypes::ACos;                 \
  using DimTypes::CosH;                 \
  using DimTypes::SinH;                 \
  using DimTypes::TanH;                 \
  using DimTypes::ACosH;                \
  using DimTypes::ASinH;                \
  using DimTypes::ATanH;                \
  using DimTypes::IsZero;               \
  using DimTypes::IsFinite;             \
  using DimTypes::IsInf;                \
  using DimTypes::IsNaN;                \
  using DimTypes::IsNeg;                \
  using DimTypes::IsPos;                \
  using DimTypes::ATan2;

