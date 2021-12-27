// Copyright 2021 E.
#ifndef UTIL_FLOAT_H_
#define UTIL_FLOAT_H_

#ifdef COMPUTE_PARTITION_MPFR
#include <boost/multiprecision/mpfr.hpp>
namespace boost {

inline void throw_exception(const std::exception& e) {
  fprintf(stderr, "Boost exception: %s\n", e.what());
  std::abort();
}

}  // namespace boost
#endif

#endif  // UTIL_FLOAT_H_
