// vim:ts=2:et
//===========================================================================//
//                            "CEMathsTest.cpp":                             //
//     Test for Our Implementation of "constexpr" Mathematical Functions     //
//===========================================================================//
#define  DIMTYPES_FORCE_OWN_ELEM_FUNCS_IMPL 1
#include "DimTypes/DimTypes.hpp"
#include <iostream>

using namespace DimTypes;

namespace
{
  //-------------------------------------------------------------------------//
  // "ErrAccum" Class:                                                       //
  //-------------------------------------------------------------------------//
  // Accumulator of Errors:
  template<typename F>
  class ErrAccum
  {
  private:
    F m_maxErr    = 0.0;
    F m_argMaxErr = NaN<F>;

  public:
    void Update(F a_x, F a_val1, F a_val0)
    {
      F base   = std::fabs(a_val0);
      F absErr = std::fabs(a_val1 - a_val0);
      // The "err" is absolute or relative, depending on the "base":
      F err    = (base >= F(1.0)) ? (absErr / base) : absErr;
      if (err > m_maxErr)
      {
        m_maxErr    = err;
        m_argMaxErr = a_x;
      }
    }

    // Get the accumulated MaxErr and ArgMaxErr:
    F GetMaxErr   () const { return m_maxErr;    }
    F GetArgMaxErr() const { return m_argMaxErr; }
  };

  //-------------------------------------------------------------------------//
  // TestFuncs:                                                              //
  //-------------------------------------------------------------------------//
  template<typename F>
  inline void TestFuncs(char const* a_title)
  {
    assert (a_title != nullptr);
    std::cout << "======= " << a_title << " =======" << std::endl;
    std::cout << "Eps="     << Eps<F>  << std::endl;

    ErrAccum<F> expErrs;
    ErrAccum<F> logErrs;
    ErrAccum<F> sqrtErrs;
    ErrAccum<F> cbrtErrs;
    ErrAccum<F> cosErrs;
    ErrAccum<F> sinErrs;
    ErrAccum<F> tanErrs;
    ErrAccum<F> atanErrs;
    ErrAccum<F> asinErrs;
    ErrAccum<F> acosErrs;
    ErrAccum<F> atan2Errs;

    for (F x = F(-80.0); x <= F(80.0); x += F(0.03125))
    {
      // Exp:
      F expX = Exp<F>(x);
      expErrs .Update(x, expX,    std::exp(x));
      
      // Log  of the previously-computed Exp:
      F logEX = Log<F>(expX);
      logErrs .Update(expX, logEX, std::log(expX));
      logErrs .Update(expX, logEX, x);

      // SqRt of the previously-computed Exp:
      F sqrtEX = SqRt<F>(expX);
      sqrtErrs.Update(expX, sqrtEX, std::sqrt(expX));
      sqrtErrs.Update(expX, sqrtEX, Exp<F>(x/F(2.0)));

      // CbRt of "x" and of the previously-computed Exp:
      F cbrtX  = CbRt<F>(x);
      cbrtErrs.Update(x,    cbrtX,  std::cbrt(x));
      F cbrtEX = CbRt<F>(expX);
      cbrtErrs.Update(expX, cbrtEX, std::cbrt(expX));
      cbrtErrs.Update(expX, cbrtEX, Exp<F>(x/F(3.0)));

      // Cos:
      F cosX  = Cos<F>(x);
      cosErrs .Update(x, cosX,  std::cos(x));

      // Sin:
      F sinX  = Sin<F>(x);
      sinErrs .Update(x, sinX,  std::sin(x));

      // Tan:
      F tanX  = Tan<F>(x);
      tanErrs .Update(x, tanX,  std::tan(x));

      // ATan:
      F atanX = ATan<F>(x);
      atanErrs.Update(x, atanX, std::atan(x));

      // Tan(ATan):
      F backX = Tan<F>(atanX);
      atanErrs.Update(atanX, backX, x);

      // ASin:
      if (Abs(x) <= F(1.0))
      {
        F asinX = ASin<F>(x);
        asinErrs.Update(x, asinX, std::asin(x));

        // Sin(ASin):
        F bsX = Sin<F>(asinX);
        asinErrs.Update(asinX, bsX, x);
      }

      // ACos:
      if (Abs(x) <= F(1.0))
      {
        F acosX = ACos<F>(x);
        acosErrs.Update(x, acosX, std::acos(x));

        // Cos(ACos):
        F bcX = Cos<F>(acosX);
        acosErrs.Update(acosX, bcX, x);
      }

      // ATan2:
      for (F y = F(-80.0); y <= F(80.0); y += F(0.03125))
      {
        F atanYX = ATan2<F>(y, x);
        if (LIKELY(std::isfinite(atanYX)))
          atan2Errs.Update(x, atanYX, std::atan2(y, x));
      }
    }
    // Results:
    std::cout << "Exp  : MaxErr=[" << expErrs  .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << expErrs  .GetArgMaxErr() << std::endl;
    std::cout << "Log  : MaxErr=[" << logErrs  .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << logErrs  .GetArgMaxErr() << std::endl;
    std::cout << "SqRt : MaxErr=[" << sqrtErrs .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << sqrtErrs .GetArgMaxErr() << std::endl;
    std::cout << "CbRt : MaxErr=[" << cbrtErrs .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << cbrtErrs .GetArgMaxErr() << std::endl;
    std::cout << "Cos  : MaxErr=[" << cosErrs  .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << cosErrs  .GetArgMaxErr() << std::endl;
    std::cout << "Sin  : MaxErr=[" << sinErrs  .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << sinErrs  .GetArgMaxErr() << std::endl;
    std::cout << "Tan  : MaxErr=[" << tanErrs  .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << tanErrs  .GetArgMaxErr() << std::endl;
    std::cout << "ATan : MaxErr=[" << atanErrs .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << atanErrs .GetArgMaxErr() << std::endl;
    std::cout << "ASin : MaxErr=[" << asinErrs .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << asinErrs .GetArgMaxErr() << std::endl;
    std::cout << "ACos : MaxErr=[" << acosErrs .GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << acosErrs .GetArgMaxErr() << std::endl;
    std::cout << "ATan2: MaxErr=[" << atan2Errs.GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="     << atan2Errs.GetArgMaxErr() << std::endl;
  }
}

int main()
{
  TestFuncs<float>      ("FLOAT      ");
  TestFuncs<double>     ("DOUBLE     ");
  TestFuncs<long double>("LONG DOUBLE");
	return 0;
}
