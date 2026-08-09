// vim:ts=2:et
//===========================================================================//
//                         "DimTypes/CompilerMacros.h":                      //
//===========================================================================//
#pragma  once

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

//===========================================================================//
// Branch Prediction:                                                        //
//===========================================================================//
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

