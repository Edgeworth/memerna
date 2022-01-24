// Copyright 2021 Eliot Courtney.
#ifndef UTIL_ERROR_H_
#define UTIL_ERROR_H_

#include <cstdio>
#include <stdexcept>

#include "util/string.h"

// Like assert, but can't be disabled.
#define verify(expr, ...)                                     \
  do {                                                        \
    if (!(expr)) {                                            \
      auto msg = ::mrna::sfmt("%s:%d: ", __func__, __LINE__); \
      msg += ::mrna::sfmt(__VA_ARGS__);                       \
      msg += '\n';                                            \
      throw ::std::runtime_error(msg);                        \
    }                                                         \
  } while (0)

#define error(...) verify(false, __VA_ARGS__)

#define bug() verify(false, "bug")

#endif  // UTIL_ERROR_H_
