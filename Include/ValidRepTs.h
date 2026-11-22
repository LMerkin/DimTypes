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
  // The Category of Types "R" which can serve as "RepT"s in "DimQ"s, implemen-
  // ted via the following predicate.
  // NB: This predicate can be further EXTENDED (unlike Concepts)  by partial
  // specialization!
  //
  template<typename R>
  constexpr inline bool IsValidDimQRepT =
    std::is_floating_point_v<R>;

  template<typename R>
  constexpr inline bool IsValidDimQRepT<std::complex<R>> =
    std::is_floating_point_v<R>;

  // Conversions of RepTypes (NB: as a Concept, there are reasons for that):
  template<typename From, typename To>
  concept DimQRepConv =
    IsValidDimQRepT<From> &&
    IsValidDimQRepT<To>   &&
    std::is_constructible_v<To, From const&>;
}
// End namespace DimTypes

