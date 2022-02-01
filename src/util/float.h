// Copyright 2021 Eliot Courtney.
#ifndef UTIL_FLOAT_H_
#define UTIL_FLOAT_H_

#ifdef USE_MPFR
#include <boost/multiprecision/mpfr.hpp>
namespace boost {

inline void throw_exception(const std::exception& e) {
  fprintf(stderr, "Boost exception: %s\n", e.what());
  std::abort();
}

}  // namespace boost
#endif  // USE_MPFR

namespace mrna {

#ifndef FLOAT_BITS
#define FLOAT_BITS 64
#endif

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
#endif  // USE_MPFR

#if FLOAT_BITS == 32
inline const flt EP{1e-3};
#elif FLOAT_BITS == 64
inline const flt EP{1e-5};  // For parition this has to be a bit smaller than normal.
#elif FLOAT_BITS == 80
inline const flt EP{1e-8};
#else
inline const flt EP{1e-30};
#endif

}  // namespace mrna

#endif  // UTIL_FLOAT_H_
