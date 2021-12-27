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
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<FLOAT_BITS>> flt;
#else
#if FLOAT_BITS == 32
typedef float flt;
#elif FLOAT_BITS == 64
typedef double flt;
#elif FLOAT_BITS == 80
typedef long double flt;
#else
static_assert(false, "unknown float bits")
#endif
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

}  // namespace mrna

#endif  // UTIL_FLOAT_H_
