// vim:ts=2:et
//===========================================================================//
//                            "CEMathsTest.cpp":                             //
//     Test for Our Implementation of "constexpr" Mathematical Functions     //
//===========================================================================//
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
    using Bits::CEMaths::Eps;

    assert (a_title != nullptr);
    std::cout << "======= " << a_title << " =======" << std::endl;
    std::cout << "Eps="     << Eps<F>  << std::endl;

    ErrAccum<F> expErrs;
    ErrAccum<F> logErrs;
    ErrAccum<F> cosErrs;
    ErrAccum<F> sinErrs;

    for (F x = F(-80.0); x <= F(80.0); x += F(0.03125))
    {
      // Exp:
      F expX = Exp(x);
      expErrs.Update(x, expX, std::exp(x));
      
      // Log of the previously-computed Exp:
      F logEX = Log(expX);
      // logErrs.Update(x, logEX, x);
      logErrs.Update(x, logEX, std::log(expX));

      // Cos:
      F cosX = Cos(x);
      cosErrs.Update(x, cosX, std::cos(x));

      // Sin:
      F sinX = Sin(x);
      sinErrs.Update(x, sinX, std::sin(x));
    }
    // Results:
    std::cout << "Exp: MaxErr=[" << expErrs.GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="   << expErrs.GetArgMaxErr() << std::endl;
    std::cout << "Log: MaxErr=[" << logErrs.GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="   << logErrs.GetArgMaxErr() << std::endl;
    std::cout << "Cos: MaxErr=[" << cosErrs.GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="   << cosErrs.GetArgMaxErr() << std::endl;
    std::cout << "Sin: MaxErr=[" << sinErrs.GetMaxErr() / Eps<F>
              << "*Eps]\t@ x="   << sinErrs.GetArgMaxErr() << std::endl;
  }
}

int main()
{
  TestFuncs<float>      ("FLOAT      ");
  TestFuncs<double>     ("DOUBLE     ");
  TestFuncs<long double>("LONG DOUBLE");
	return 0;
}
