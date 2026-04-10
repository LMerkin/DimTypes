// vim:ts=2:et
//===========================================================================//
//                          "DimTypes/DimTypes.hpp":                         //
//          Types for Physical Dimensions and Dimensioned Quantities         //
//===========================================================================//
#pragma  once
#include "Bits/Encodings.hpp"
#include "Bits/Macros.h"

namespace DimTypes
{
  //=========================================================================//
  // Lifted "constexpr" templated mathematical Consts and Functions:         //
  //=========================================================================//
  using namespace Bits::CEMaths;

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
    // Non-Default Ctor: Lifting a value into "DimQ" (by default, 0). NB:  We
    // DELIBERATELY DO NOT INFER the "RepT" from the type of the Ctor arg, as
    // this may cause ambiguities in some cases:
    //
    constexpr          DimQ()            : m_val(0.0)         {}

    // "R" is required to be a floating-point type to avoid a conflict with the
    // Copy Ctor below (where "R" is "DimQ" itself):
    template<typename  R>
    constexpr explicit DimQ(R const& val) requires(std::is_floating_point_v<R>)
    : m_val(RepT(val)) {}

    // Copy Ctor:
    constexpr DimQ(DimQ const& a_right): m_val(a_right.m_val) {}

    // Assignment:
    constexpr DimQ& operator= (DimQ const& a_right)
    {
      m_val = a_right.m_val;
      return *this;
    }

    // If the assignment LHS is DimLess, the RHS may be of any floating-point
    // type "R"   (again, this constraint is required to avoid confusion with
    // the "DimQ" itself):
    template<typename R>
    constexpr DimQ operator= (R const& a_right)
    requires(E == 0 && std::is_floating_point_v<R>)
    {
      m_val = RepT(a_right);
      return *this;
    }

    // Conversion FROM another Rep (but E, U and MaxDims must match);
    // this is an "Extended Copy Ctor":
    template<typename  FromR>
    constexpr explicit DimQ(DimQ<E, U, FromR, MaxDims> const& a_right)
      : m_val(RepT(a_right.Magnitude()))
      {}

    // Conversion TO another Rep (again, E, U and MaxDims are unchanged);
    // this is just a shortcut to the above "Extended Copy Ctor":
    template<typename ToR>
    constexpr friend DimQ<E, U, ToR, MaxDims> ToRepT(DimQ const& a_right)
      { return DimQ<E, U, ToR, MaxDims>(a_right); }

    // The dimensioned unit for the given dimensioned quantity. NB: This method
    // is not a ctor: it takes an existing "DimQ" and returns a similar one but
    // with the unitary magnitude:
    //
    constexpr DimQ UnitOf() const           { return DimQ(RepT(1.0));  }

    // Getting the exponent code (this is primarily for testing):
    constexpr uint64_t GetDimsCode () const { return E; }

    // Getting the Units code (also for testing):
    constexpr uint64_t GetUnitsCode() const { return U; }

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

    // NB:  In addition, for some types (eg Angles expressed in Radians), the
    // user may allow direct conversion into "RepT" by specialising the above
    // operator method...

    //=======================================================================//
    // Arithmetic Operations -- Dimensions Unchanged:                        //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Addition and Subtraction of "DimQ"s:                                  //
    //-----------------------------------------------------------------------//
    // This is only possible if the dimensions and RepT are same, and the units
    // can be unified (i.e. same units for non-0 exponents),  in which case the
    // 1st arg's units can still be used:
    //
    // Addition:
    template<uint64_t F, uint64_t V>
    constexpr DimQ operator+   (DimQ<F,V,RepT,MaxDims> const& a_right) const
    {
      static_assert(E == F,             "ERROR: +: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: +: Units do not unify");
      return DimQ(m_val + a_right.m_val);
    }

    template<uint64_t F, uint64_t V>
    constexpr DimQ& operator+= (DimQ<F,V,RepT,MaxDims> const& a_right)
    {
      static_assert(E == F,             "ERROR: +=: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: +=: Units do not unify");
      m_val += a_right.m_val;
      return *this;
    }

    // Subtraction of "DimQs":
    // Same constraints as for addition:
    template<uint64_t F, uint64_t V>
    constexpr DimQ  operator- (DimQ<F,V,RepT,MaxDims> const& a_right) const
    {
      static_assert(E == F,             "ERROR: -: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: -: Units do not unify");
      return DimQ(m_val - a_right.m_val);
    }

    template<uint64_t F, uint64_t V>
    constexpr DimQ& operator-=(DimQ<F,V,RepT,MaxDims> const& a_right)
    {
      static_assert(E == F,             "ERROR: -=: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: -=: Units do not unify");
      m_val -= a_right.m_val;
      return *this;
    }

    //-----------------------------------------------------------------------//
    // Multiplication and Division by a "RepT":                              //
    //-----------------------------------------------------------------------//
    // NB: the following method can be used to create new dimensioned values
    // from the fundamental units. XXX: when used this way, it is slighly in-
    // efficient because of multiplication by 1 and extra copying;   but the
    // impact of that should be negligible:
    //
    // Multiplication:
    constexpr  DimQ operator*   (RepT const& a_right) const
      { return DimQ(m_val * a_right); }

    constexpr  DimQ& operator*= (RepT const& a_right)
    {
      m_val *= a_right;
      return *this;
    }

    constexpr friend DimQ operator* (RepT const& a_left, DimQ const& a_right)
      { return DimQ(a_left * a_right.m_val); }

    // In particular, Unary Negation:
    constexpr  DimQ operator- () const
      { return DimQ(- m_val); }

    // Division by a "RepT":
    constexpr  DimQ operator/   (RepT const& a_right) const
      { return DimQ(m_val / a_right); }

    constexpr  DimQ& operator/= (RepT const& a_right)
    {
      m_val /= a_right;
      return *this;
    }

    //-----------------------------------------------------------------------//
    // "Abs", "Floor", "Ceil", "Round": "DimQ" -> "DimQ":                    //
    //-----------------------------------------------------------------------//
    // XXX: Currently for Real args only (even "Abs"):
    //
#   ifdef  DIMQ_ROUND_FUNC
#   undef  DIMQ_ROUND_FUNC
#   endif
#   define DIMQ_ROUND_FUNC(FuncName) \
    constexpr  friend DimQ FuncName(DimQ const& a_right) \
    requires(!IsComplex<RepT>) \
      { return DimQ(FuncName(a_right.m_val)); }

    DIMQ_ROUND_FUNC(Abs)
    DIMQ_ROUND_FUNC(Floor)
    DIMQ_ROUND_FUNC(Ceil)
    DIMQ_ROUND_FUNC(Round)
#   undef  DIMQ_ROUND_FUNC

    //=======================================================================//
    // Arithmetic operations which result in Dimensions change:              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Multiplication of "DimQ"s:                                            //
    //-----------------------------------------------------------------------//
    // Dimension exponents are added (or subtracted) up; units must unify;
    // afterwards, units of zero-exp dims are reset, otherwise false type-che-
    // cking failures may occur:
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

    constexpr friend DimQ<En::SubExp(0UL,E), U, RepT, MaxDims>
    operator/ (RepT const& a_left, DimQ const& a_right)
    {
      return DimQ<En::SubExp(0UL,E), U, RepT, MaxDims>
                 (a_left / a_right.m_val);
    }

    //-----------------------------------------------------------------------//
    // Integral and Rational Powers:                                         //
    //-----------------------------------------------------------------------//
    // the power must be given by constant exprs
    // so it is analysable at compile time:
    // NB: for applying these methods from outside "DimQ",  explicit  "template
    // pow" expr is required, otherwise the template is not syntactically recog-
    // nised!
    // NB: "CleanUpUnits" is required in case if M==0:
    //-----------------------------------------------------------------------//
    // "IPow": Integral Power, as Prefix Func:                               //
    //-----------------------------------------------------------------------//
    template<int M>
    constexpr friend
    DimQ<En::MultExp(E,M), En::CleanUpUnits(En::MultExp(E,M), U),
         RepT,             MaxDims>
    IPow(DimQ const& a_right)
    {
      return
        DimQ<En::MultExp(E,M), En::CleanUpUnits(En::MultExp(E,M), U),
             RepT,             MaxDims>
            (En::template IntPow<M>(a_right.m_val));
    }

    //-----------------------------------------------------------------------//
    // "Sqr" and "Cube" are particularly-important cases of "IPow":          //
    //-----------------------------------------------------------------------//
    // Provided as Prefix Funcs:
    //
    constexpr friend
    DimQ<En::MultExp(E,2), En::CleanUpUnits(En::MultExp(E,2), U), RepT, MaxDims>
    Sqr(DimQ const& a_right)
      { return IPow<2>(a_right); }

    constexpr friend
    DimQ<En::MultExp(E,3), En::CleanUpUnits(En::MultExp(E,3), U), RepT, MaxDims>
    Cube(DimQ const& a_right)
      { return IPow<3>(a_right); }

    //-----------------------------------------------------------------------//
    // "RPow": (General) Rational Power:                                     //
    //-----------------------------------------------------------------------//
    // XXX: In general it is NOT "constexpr" for C++ < 26, but it is "constexpr"
    // if it can be reduced to a composition of "SqRt" and "CbRt":
    //
    template<int M, int N>
    constexpr friend
    DimQ<En::DivExp(En::MultExp(E,M),N),
         En::CleanUpUnits(En::DivExp(En::MultExp(E,M),N), U),
         RepT,
         MaxDims>
    RPow(DimQ const& a_right)
    {
      // Multiples of "PMod" (including 0) are uninvertible in the curr rep;
      // and for a fractional power, "constexpr" is not possible yet:
      //
      static_assert(N > 0 && N % En::IPMod != 0, "RPow: Invalid Denom");
      return
        DimQ<En::DivExp(En::MultExp(E,M),N),
             En::CleanUpUnits(En::DivExp(En::MultExp(E,M),N), U),
             RepT,
             MaxDims>
            (En::template FracPow<M, N>(a_right.m_val));
    }

    //-----------------------------------------------------------------------//
    // Shortcuts: "SqRt" and "CbRt", as Prefix Funcs:                        //
    //-----------------------------------------------------------------------//
    // Here M==1, so "CleanUpUnits" is not required:
    //
    constexpr friend DimQ<En::DivExp(E,2), U, RepT, MaxDims> SqRt
      (DimQ const& a_right)
      { return     RPow<1,2>(a_right); }

    constexpr friend DimQ<En::DivExp(E,3), U, RepT, MaxDims> CbRt
      (DimQ const& a_right)
      { return     RPow<1,3>(a_right); }

    //-----------------------------------------------------------------------//
    // "ATan2" on "DimQs": Prefix Func:                                      //
    //-----------------------------------------------------------------------//
    template<uint64_t F,  uint64_t V>
    constexpr friend RepT ATan2
      (DimQ const& a_y, DimQ<F,V,RepT,MaxDims> const& a_x)
    requires(!IsComplex<RepT>)
    {
      static_assert(E == F,             "ERROR: ATan2: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: ATan2: Units do not unify");
      return ATan2(a_y.m_val, a_x.m_val);
    }

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
#   define DIMQ_CMP(Op)  \
    template<uint64_t V> \
    constexpr bool operator Op(DimQ<E, V, RepT, MaxDims> const& a_right) const \
    { \
      static_assert(En::UnitsOK(E,U,V), "ERROR: Units do not unify"); \
      return m_val Op a_right.m_val; \
    }
    DIMQ_CMP(==)
    DIMQ_CMP(!=)
    DIMQ_CMP(>)
    DIMQ_CMP(>=)
    DIMQ_CMP(<)
    DIMQ_CMP(<=)
#   undef  DIMQ_CMP

    //-----------------------------------------------------------------------//
    // Approximate Equality:                                                 //
    //-----------------------------------------------------------------------//
    // "ApproxEquals" is a method:
    template<uint64_t V>
    constexpr bool ApproxEquals
    (
      DimQ<E, V, RepT, MaxDims> const& a_right,
      RepT                      const& a_tol = DefaultTol<RepT>
    )
    const
    {
      static_assert(En::UnitsOK(E,U,V), "ERROR: Units do not unify");
      return ApproxEqual(m_val, a_right.m_val, a_tol);
    }

    // "ApproxEqual" (w/o "s"!) is a friend (similar to "ApproxEqual" on "RepT",
    // so the same call API can be used for both):
    template<uint64_t V>
    constexpr friend bool ApproxEqual
    (
      DimQ<E, U, RepT, MaxDims> const& a_left,
      DimQ<E, V, RepT, MaxDims> const& a_right,
      RepT                      const& a_tol = DefaultTol<RepT>
    )
    { return a_left.ApproxEquals(a_right, a_tol); }
  };
  // End "DimQ" class

  //=========================================================================//
  // NameSpace-Level Functions:                                              //
  //=========================================================================//
  // IMPORTANT: DO NOT make them "DimQ friends", as this would limit the templ-
  // ate substitution / overloading resolution:
  //-------------------------------------------------------------------------//
  // Unary Mathematical Funcs: DimLess->RepT:                                //
  //-------------------------------------------------------------------------//
# ifdef  DIMLESS_UNARY_FUNC
# undef  DIMLESS_UNARY_FUNC
# endif
# define DIMLESS_UNARY_FUNC(FuncName) \
  template<uint64_t E, uint64_t U, typename RepT, unsigned MaxDims> \
  constexpr RepT FuncName(DimQ<E, U, RepT, MaxDims> const& a_right) \
  { \
    static_assert(E == 0, "ERROR: Must be DimLess"); \
    return FuncName(a_right.Magnitude()); \
  }
  DIMLESS_UNARY_FUNC(Exp)
  DIMLESS_UNARY_FUNC(Log)
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
  template<uint64_t E, uint64_t U, typename RepT, unsigned MaxDims> \
  constexpr bool Pred(DimQ<E, U, RepT, MaxDims> const& a_right) \
  requires(!IsComplex<RepT>) \
    { return Pred(a_right.Magnitude());  }

  DIMQ_PRED(IsZero)
  DIMQ_PRED(IsFinite)
  DIMQ_PRED(IsInf)
  DIMQ_PRED(IsNaN)
  DIMQ_PRED(IsNeg)
  DIMQ_PRED(IsPos)
# undef  DIMQ_PRED
}
// End namespace DimTypes
