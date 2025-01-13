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
    RepT m_val;  // The actual value (magnitude)

    // Encodings to be used:
    using En =   Bits::Encodings<RepT, MaxDims>;

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

    // Conversion from another Rep (but E and U must match):
    template<typename R>
    constexpr explicit DimQ(DimQ<E,U,R,MaxDims> const& a_right)
      : m_val(RepT(a_right.m_val))
      {}

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
    constexpr DimQ operator+   (DimQ<F,V,RepT,MaxDims> a_right) const
    {
      static_assert(E == F,             "ERROR: +: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: +: Units do not unify");
      return DimQ(m_val + a_right.m_val);
    }

    template<uint64_t F, uint64_t V>
    constexpr DimQ& operator+= (DimQ<F,V,RepT,MaxDims> a_right)
    {
      static_assert(E == F,             "ERROR: +=: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: +=: Units do not unify");
      m_val += a_right.m_val;
      return *this;
    }

    // Subtraction of "DimQs":
    // Same constraints as for addition:
    template<uint64_t F, uint64_t V>
    constexpr DimQ operator- (DimQ<F,V,RepT,MaxDims>  a_right) const
    {
      static_assert(E == F,             "ERROR: -: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: -: Units do not unify");
      return DimQ(m_val - a_right.m_val);
    }

    template<uint64_t F, uint64_t V>
    constexpr DimQ& operator-=(DimQ<F,V,RepT,MaxDims> a_right)
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
    constexpr  DimQ operator*(RepT a_right) const
      { return DimQ(m_val * a_right); }

    constexpr DimQ& operator*= (RepT a_right)
    {
      m_val *= a_right;
      return *this;
    }

    constexpr friend DimQ operator* (RepT a_left, DimQ a_right)
      { return DimQ(a_left * a_right.m_val); }

    // In particular, Unary Negation:
    constexpr  DimQ operator- () const
      { return DimQ(- m_val); }

    // Division by a "RepT":
    constexpr  DimQ operator/   (RepT a_right) const
      { return DimQ(m_val / a_right); }

    constexpr  DimQ& operator/= (RepT a_right)
    {
      m_val /= a_right;
      return *this;
    }

    //-----------------------------------------------------------------------//
    // "Abs", "Floor", "Ceil", "Round":                                      //
    //-----------------------------------------------------------------------//
    // XXX: They are currently provided for Real "RepT" only (for Complex
    // "RepT", the corresp ops will cause compilation errors):
    //
    constexpr  DimQ Abs  () const { return DimQ(Bits::CEMaths::Abs  (m_val)); }
    constexpr  DimQ Floor() const { return DimQ(Bits::CEMaths::Floor(m_val)); }
    constexpr  DimQ Ceil () const { return DimQ(Bits::CEMaths::Ceil (m_val)); }
    constexpr  DimQ Round() const { return DimQ(Bits::CEMaths::Round(m_val)); }

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
    operator* (DimQ<F, V, RepT, MaxDims> a_right) const
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
    operator/ (DimQ<F, V, RepT, MaxDims> a_right) const
    {
      return  DimQ<En::SubExp(E,F),
                   En::CleanUpUnits(En::SubExp(E,F), En::UnifyUnits(E,F,U,V)),
                   RepT,
                   MaxDims>
             (m_val / a_right.Magnitude());
    }

    constexpr friend DimQ<En::SubExp(0UL,E), U, RepT, MaxDims>
    operator/ (RepT a_left, DimQ a_right)
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
    // "IPow": Integral Power:                                               //
    //-----------------------------------------------------------------------//
    template<int M>
    constexpr DimQ<En::MultExp(E,M), En::CleanUpUnits(En::MultExp(E,M), U),
                   RepT,             MaxDims>
    IPow() const
    {
      return
        DimQ<En::MultExp(E,M), En::CleanUpUnits(En::MultExp(E,M), U),
             RepT,             MaxDims>
            (En::template IntPow<M>(m_val));
    }

    //-----------------------------------------------------------------------//
    // "Sqr" and "Cube" are particularly-important case of "IPow":           //
    //-----------------------------------------------------------------------//
    constexpr DimQ<En::MultExp(E,2), En::CleanUpUnits(En::MultExp(E,2), U),
                   RepT,             MaxDims>
    Sqr()  const
      { return IPow<2>(); }

    constexpr DimQ<En::MultExp(E,3), En::CleanUpUnits(En::MultExp(E,3), U),
                   RepT,             MaxDims>
    Cube() const
      { return IPow<3>(); }

    //-----------------------------------------------------------------------//
    // "RPow": (General) Rational Power:                                     //
    //-----------------------------------------------------------------------//
    // XXX: In general it is NOT "constexpr" for C++ < 26, but it is "constexpr"
    // if it can be reduced to a composition of "SqRt" and "CbRt":
    //
    template<int M, int N>
    constexpr DimQ<En::DivExp(En::MultExp(E,M),N),
                   En::CleanUpUnits(En::DivExp(En::MultExp(E,M),N), U),
                   RepT,
                   MaxDims>
    RPow() const
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
            (En::template FracPow<M, N>(m_val));
    }

    //-----------------------------------------------------------------------//
    // Shortcuts: "SqRt" and "CbRt":                                         //
    //-----------------------------------------------------------------------//
    // Here M==1, so "CleanUpUnits" is not required:
    constexpr DimQ<En::DivExp(E,2), U, RepT, MaxDims> SqRt() const
      { return RPow<1,2>(); }

    constexpr DimQ<En::DivExp(E,3), U, RepT, MaxDims> CbRt() const
      { return RPow<1,3>(); }

    //-----------------------------------------------------------------------//
    // Unary DimLess Elementary Functions:                                   //
    //-----------------------------------------------------------------------//
#   ifdef  DIMLESS_UNARY_FUNC
#   undef  DIMLESS_UNARY_FUNC
#   endif
#   define DIMLESS_UNARY_FUNC(FuncName) \
    constexpr DimQ<0, 0, RepT, MaxDims> FuncName(DimQ a_x) const \
    { \
      static_assert(E==0, "ERROR: " #FuncName ": Must be DimLess"); \
      return  DimQ<0, 0, RepT, MaxDims>(Bits::CEMaths::FuncName(a_x.m_val)); \
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
#   undef DIMLESS_UNARY_FUNC

    //-----------------------------------------------------------------------//
    // "ATan2" on "DimQs":                                                   //
    //-----------------------------------------------------------------------//
    // "y" is *this, "x" is the explicit arg:
    //
    template<uint64_t F, uint64_t V>
    constexpr RepT ATan2(DimQ<F,V,RepT,MaxDims> a_x) const
    {
      static_assert(E == F,             "ERROR: ATan2: Different Dims");
      static_assert(En::UnitsOK(E,U,V), "ERROR: ATan2: Units do not unify");
      return Bits::CEMaths::ATan2(m_val, a_x.m_val);
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
    constexpr bool operator Op (DimQ<E, V, RepT, MaxDims> a_right) const \
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

    constexpr bool IsZero  () const  { return m_val == RepT(0.0);   }
    constexpr bool IsFinite() const  { return std::isfinite(m_val); }
    constexpr bool IsNaN   () const  { return std::isnan   (m_val); }

    // NB: The following methods will not compile if the field "RepT" is not
    // ordered (eg for complex numbers). However, if they  are not  actually
    // used, this will cause no harm:
    //
    constexpr bool IsNeg   () const  { return m_val <  RepT(0.0); }
    constexpr bool IsPos   () const  { return m_val >  RepT(0.0); }

    // Approximate Equality:
    template<uint64_t V>
    constexpr bool ApproxEquals
    (
      DimQ<E, V, RepT, MaxDims> a_right,
      RepT                      a_tol = Bits::CEMaths::DefaultTol<RepT>
    )
    const
    {
      static_assert(En::UnitsOK(E,U,V), "ERROR: Units do not unify");
      return Bits::CEMaths::ApproxEqual(m_val, a_right.m_val, a_tol);
    }
  };
  // End "DimQ" class

  //=========================================================================//
  // Syntactic Sugar: Above methods in Prefix notation:                      //
  //=========================================================================//
  // XXX: INCREDIBLY, in C++20,  these functions are visible directly  (ie w/o
  // the "DimTypes::" prefix!!!) outside the "DimTypes" namespace, if the cor-
  // resp arg "DimQ" obj is there.  So their visibility is very much like that
  // of "DimQ" methods:
  //
  template <uint64_t E, uint64_t U, typename RepT, unsigned   MaxDims>
  constexpr DimQ<E, U, RepT, MaxDims> UnitOf(DimQ<E, U, RepT, MaxDims> a_dimq)
    { return a_dimq.UnitOf(); }

  template <uint64_t E, uint64_t U, typename RepT, unsigned  MaxDims>
  constexpr uint64_t GetDimsCode (DimQ<E, U, RepT, MaxDims>  a_dimq)
    { return a_dimq .GetDimsCode(); }

  template <uint64_t E, uint64_t U, typename RepT, unsigned  MaxDims>
  constexpr uint64_t GetUnitsCode(DimQ<E, U, RepT, MaxDims>  a_dimq)
    { return a_dimq.GetUnitsCode(); }

  template <uint64_t E, uint64_t U, typename RepT, unsigned  MaxDims>
  constexpr RepT     Magnitude   (DimQ<E, U, RepT, MaxDims>  a_dimq)
    { return a_dimq .Magnitude();   }

  template <uint64_t E, uint64_t U, typename RepT, unsigned  MaxDims>
  constexpr DimQ<E, U, RepT, MaxDims> Abs  (DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.Abs();   }

  template< uint64_t E, uint64_t U, typename RepT, unsigned  MaxDims>
  constexpr DimQ<E, U, RepT, MaxDims> Floor(DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.Floor(); }

  template< uint64_t E, uint64_t U, typename RepT, unsigned  MaxDims>
  constexpr DimQ<E, U, RepT, MaxDims> Ceil (DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.Ceil();  }

  template< uint64_t E,  uint64_t U, typename RepT, unsigned MaxDims>
  constexpr DimQ<E, U, RepT, MaxDims> Round(DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.Round(); }

  template<int M,  uint64_t E,  uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto IPow(DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.template IPow<M>(); }

  template< uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto Sqr(DimQ<E, U, RepT, MaxDims>  a_right)
    { return a_right.Sqr(); }

  template< uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto Cube(DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.Cube(); }

  template<int M,  int N, uint64_t E, uint64_t U,
           typename RepT, unsigned MaxDims>
  constexpr auto RPow(DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.template RPow<M,N>(); }

  template< uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto SqRt(DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.SqRt(); }

  template< uint64_t E, uint64_t U, typename RepT, unsigned MaxDims>
  constexpr auto CbRt(DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.CbRt(); }

  template< uint64_t E, uint64_t  U, typename RepT, unsigned MaxDims>
  constexpr bool IsZero  (DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.IsZero();    }

  template< uint64_t E, uint64_t  U, typename RepT, unsigned MaxDims>
  constexpr bool IsFinite(DimQ<E, U,RepT, MaxDims>  a_right)
    { return a_right.IsFinite();  }

  template< uint64_t E, uint64_t  U, typename RepT, unsigned MaxDims>
  constexpr bool IsNeg   (DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.IsNeg();     }

  template< uint64_t E, uint64_t  U, typename RepT, unsigned MaxDims>
  constexpr bool IsPos   (DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.IsPos();     }

  template< uint64_t E, uint64_t  U, typename RepT, unsigned MaxDims>
  constexpr bool IsNaN   (DimQ<E, U, RepT, MaxDims> a_right)
    { return a_right.IsNaN();     }

  template< uint64_t E, uint64_t  U, typename RepT, unsigned MaxDims>
  constexpr RepT ATan2
    (DimQ<E, U, RepT, MaxDims> a_y, DimQ<E, U, RepT, MaxDims> a_x)
    { return a_y.ATan2(a_x); }

  //=========================================================================//
  // Lifted "constexpr" Mathematical Functions:                              //
  //=========================================================================//
  using Bits::CEMaths::NaN;
  using Bits::CEMaths::Inf;
  using Bits::CEMaths::Eps;
  using Bits::CEMaths::DefaultTol;

  using Bits::CEMaths::SqRt2;
  using Bits::CEMaths::SqRt1_2;
  using Bits::CEMaths::SqRt3;

  using Bits::CEMaths::Pi;
  using Bits::CEMaths::TwoPi;
  using Bits::CEMaths::Pi_2;
  using Bits::CEMaths::Pi_4;

  using Bits::CEMaths::Abs;
  using Bits::CEMaths::Floor;
  using Bits::CEMaths::Ceil;
  using Bits::CEMaths::Round;
  using Bits::CEMaths::ApproxEqual;

  using Bits::CEMaths::Sqr;
  using Bits::CEMaths::Cube;
  using Bits::CEMaths::SqRt;
  using Bits::CEMaths::CbRt;
  using Bits::CEMaths::Pow;

  using Bits::CEMaths::CosSin;
  using Bits::CEMaths::ATan2;

  //-------------------------------------------------------------------------//
  // The following functions are available on both "RepT" and DimLess qtys:  //
  //-------------------------------------------------------------------------//
# ifdef  DIMLESS_UNARY_FUNC
# undef  DIMLESS_UNARY_FUNC
# endif
# define DIMLESS_UNARY_FUNC(FuncName) \
  using Bits::CEMaths::FuncName;      \
  template<uint64_t U, typename RepT, unsigned MaxDims> \
  DimQ<0, 0, RepT, MaxDims> FuncName(DimQ<0, U, RepT, MaxDims> a_x) \
    { return a_x.FuncName(); }

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
}
// End namespace DimTypes
