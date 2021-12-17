// Copyright 2016 E.
#ifndef COMMON_H_
#define COMMON_H_

#include <cassert>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#ifdef COMPUTE_PARTITION_MPFR
#include <boost/multiprecision/mpfr.hpp>
namespace boost {

inline void throw_exception(const std::exception& e) {
  fprintf(stderr, "Boost exception: %s\n", e.what());
  std::abort();
}

}  // namespace boost
#endif

// Like assert, but can't be disabled.
#define verify(expr, ...)                             \
  do {                                                \
    if (!(expr)) {                                    \
      fprintf(stderr, "%s:%d: ", __func__, __LINE__); \
      fprintf(stderr, __VA_ARGS__);                   \
      fprintf(stderr, "\n");                          \
      exit(1);                                        \
    }                                                 \
  } while (0)

#endif  // COMMON_H_
