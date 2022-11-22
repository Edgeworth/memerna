// Copyright 2021 Eliot Courtney.
#ifndef UTIL_FLOAT_H_
#define UTIL_FLOAT_H_

#include <algorithm>
#include <cmath>
#include <cstdint>

#ifdef USE_MPFR
#include <boost/multiprecision/mpfr.hpp>
namespace boost {

inline void throw_exception(const std::exception& e) {
  std::cerr << "Boost exception: " << e.what() << '\n';
  std::abort();
}

}  // namespace boost
#endif  // USE_MPFR

namespace mrna {

constexpr int powi(int base, int exp) {
  int result = 1;
  for (int i = 0; i < exp; ++i) result *= base;
  return result;
}

#ifdef USE_MPFR
using flt = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<FLOAT_BITS>>;
#else
#if FLOAT_BITS == 32
using flt = float;
#elif FLOAT_BITS == 64
using flt = double;
#elif FLOAT_BITS == 80
using flt = long double;
#else
static_assert(false, "unknown float bits")
#endif

inline flt fabs(flt v) { return std::abs(v); }

#endif  // USE_MPFR

#if FLOAT_BITS == 32
inline const flt EP{1e-3};
#elif FLOAT_BITS == 64
inline const flt EP{1e-6};
#elif FLOAT_BITS == 80
inline const flt EP{1e-8};
#else
inline const flt EP{1e-30};
#endif

inline bool abs_eq(flt a, flt b, flt ep = EP) { return fabs(a - b) < ep; }

inline bool rel_eq(flt a, flt b, flt rel = EP) {
  return fabs(a - b) <= rel * std::max(fabs(a), fabs(b));
}

}  // namespace mrna

#endif  // UTIL_FLOAT_H_
