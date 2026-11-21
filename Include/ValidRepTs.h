// vim:ts=2:et
//===========================================================================//
//                       "DimTypes/ValidRepTs.h":                            //
//===========================================================================//
#pragma  once
#include <type_traits>
#include <complex>

namespace DimTypes
{
  //=========================================================================//
  // Representation Types:                                                   //
  //=========================================================================//
  // The Category of Types which can serve as "RepT"s in "DimQ"s (below), impl-
  // elemnted via the following predicate:
  //
  template<typename R>
  constexpr inline bool IsValidDimQRepT =
    std::is_floating_point_v<R>;

  template<typename R>
  constexpr inline bool IsValidDimQRepT<std::complex<R>> =
    std::is_floating_point_v<R>;

  // NB: This predicate can be further EXTENDED (unlike Concepts) by partial
  // specialization!
}
// End namespace DimTypes

