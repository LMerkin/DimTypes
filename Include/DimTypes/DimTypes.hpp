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
  template<uint64_t E, uint64_t U, typename RepT = double>
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
    template<uint64_t V>
    constexpr DimQ operator+   (DimQ<E,V> a_right) const
    {
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify");
      return DimQ(m_val + a_right.Magnitude());
    }
  
    template<uint64_t V>
    constexpr DimQ& operator+= (DimQ<E,V> a_right)
    {
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify");
      m_val += a_right.Magnitude();
      return *this;
    }
  
    // Subtraction of "DimQs":
    // Same constraints as for addition:
    template<uint64_t V>
    constexpr DimQ operator- (DimQ<E,V,RepT> a_right) const
    {
      static_assert(UnitsOK(E,U,V), "ERROR: Units do not unify");
      return DimQ(m_val - a_right.Magnitude());
    }
  
    template<uint64_t V>
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
    // "Abs":
    //-----------------------------------------------------------------------//
    // XXX: It is currently provided for Real "RepT" only (for Complex "RepT",
    // comparison with 0 will cause a compilation error), because Complex sem-
    // antics is quite different anyway:
    //
    constexpr  DimQ Abs  () const
      { return DimQ((m_val < RepT(0.0)) ? (-m_val): m_val); }

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
    template<uint64_t F, uint64_t V>
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
    template<uint64_t F, uint64_t V>
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
            (IntPow<RepT, M>(m_val));
    }

    // "Sqr" and "Cube" are particularly-important case of "IPow":
    constexpr DimQ<MultExp(E,2), CleanUpUnits(MultExp(E,2), U), RepT>
    Sqr()  const
      { return IPow<2>(); }

    constexpr DimQ<MultExp(E,3), CleanUpUnits(MultExp(E,3), U), RepT>
    Cube() const
      { return IPow<3>(); }

    // "RPow": (General) Rational Power, NOT "constexpr" for C++ < 26:
    template<int M, int N>
    constexpr DimQ<DivExp(MultExp(E,M),N),
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
            (FracPow<M,N,RepT>(m_val));
    }

    // Shortcuts: "SqRt" and "CbRt". Here M==1, so "CleanUpUnits" is not
    // required; again, these functions are NOT "constexpr" for C++ < 26:
    //
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
#   define DIMQ_CMP(Op)  \
    template<uint64_t V> \
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
    // ordered (eg for complex numbers). However, if they  are not  actually
    // used, this will cause no harm:
    //
    constexpr bool IsNeg   () const  { return m_val <  RepT(0.0); }
    constexpr bool IsPos   () const  { return m_val >  RepT(0.0); }
  };

  //=========================================================================//
  // Syntactic Sugar: Above methods in Prefix notation:                      //
  //=========================================================================//
  // XXX: INCREDIBLY, in C++20,  these functions are visible directly  (ie w/o
  // the "DimTypes::" prefix!!!) outside the "DimTypes" namespace, if the cor-
  // resp arg "DimQ" obj is there.  So their visibility is very much like that
  // of "DimQ" methods:
  //
  template<uint64_t E, uint64_t U, typename RepT>
  constexpr DimQ<E,U,RepT>  UnitOf(DimQ<E,U,RepT> a_dimq)
    { return a_dimq.UnitOf(); }

  template<uint64_t E, uint64_t U, typename RepT>
  constexpr uint64_t GetDimsCode  (DimQ<E,U,RepT> a_dimq)
    { return a_dimq.GetDimsCode();  }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr uint64_t GetUnitsCode(DimQ<E,U,RepT> a_dimq)
    { return a_dimq.GetUnitsCode(); }
  
  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr RepT Magnitude(DimQ<E,U,RepT> a_dimq)
    { return a_dimq.Magnitude();    }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<E,U,RepT>  Abs  (DimQ<E,U,RepT> a_right)
    { return a_right.Abs();    }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<E,U,RepT>  Floor(DimQ<E,U,RepT> a_right)
    { return a_right.Floor();  }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<E,U,RepT>  Ceil (DimQ<E,U,RepT> a_right)
    { return a_right.Ceil();  }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<E,U,RepT>  Round(DimQ<E,U,RepT> a_right)
    { return a_right.Round(); }

  template<int M,  uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<MultExp(E,M), CleanUpUnits(MultExp(E,M),U), RepT>
  IPow(DimQ<E,U,RepT> a_right)
    { return a_right.template IPow<M>(); }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<MultExp(E,2), CleanUpUnits(MultExp(E,2),U), RepT>
  Sqr(DimQ<E,U,RepT> a_right)
    { return a_right.Sqr(); }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<MultExp(E,3), CleanUpUnits(MultExp(E,3),U), RepT>
  Cube(DimQ<E,U,RepT> a_right)
    { return a_right.Cube(); }

  template<int M, int N,  uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<DivExp(MultExp(E,M),N),
                 CleanUpUnits(DivExp(MultExp(E,M),N), U), RepT>
  RPow(DimQ<E,U,RepT> a_right)
    { return a_right.template RPow<M,N>(); }

  // Again, "SqRt" and "CbRt" are NOT "constexpr" for C++ < 26:
  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<DivExp(E,2),U,RepT> SqRt(DimQ<E,U,RepT> a_right)
    { return a_right.SqRt(); }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr DimQ<DivExp(E,3),U,RepT> CbRt(DimQ<E,U,RepT> a_right)
    { return a_right.CbRt(); }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr bool IsZero  (DimQ<E,U,RepT> a_right)
    { return a_right.IsZero();     }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr bool IsFinite(DimQ<E,U,RepT> a_right)
    { return a_right.IsFinite();   }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr bool IsNeg   (DimQ<E,U,RepT> a_right)
    { return a_right.IsNeg(); }

  template< uint64_t E,  uint64_t U, typename RepT>
  constexpr bool IsPos   (DimQ<E,U,RepT> a_right)
    { return a_right.IsPos(); }

  //=========================================================================//
  // Lifetd "constexpr" Mathematical Functions:                              //
  //=========================================================================//
  using Bits::CEMaths::NaN;
  using Bits::CEMaths::Pi;

  using Bits::CEMaths::Exp;
  using Bits::CEMaths::Log;
  using Bits::CEMaths::Cos;
  using Bits::CEMaths::Sin;
  using Bits::CEMaths::CosSin;
  using Bits::CEMaths::SqRt;
  using Bits::CEMaths::CbRt;
  using Bits::CEMaths::Pow;
}
// End namespace DimTypes
