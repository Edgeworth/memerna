// Copyright 2021 E.
#ifndef UTIL_ERROR_H_
#define UTIL_ERROR_H_

#include <cstdio>
#include <stdexcept>

#include "util/string.h"

#ifdef USE_BOOST
#include <boost/stacktrace.hpp>
// Like assert, but can't be disabled.
#define verify(expr, ...)                                                   \
  do {                                                                      \
    if (!(expr)) {                                                          \
      auto msg = ::mrna::sfmt("%s:%d: ", __func__, __LINE__);               \
      msg += ::mrna::sfmt(__VA_ARGS__);                                     \
      msg += '\n';                                                          \
      msg += boost::stacktrace::to_string(boost::stacktrace::stacktrace()); \
      throw ::std::runtime_error(msg);                                      \
    }                                                                       \
  } while (0)
#else
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
#endif

#define error(...) verify(false, __VA_ARGS__)

#define bug() verify(false, "bug")

#endif  // UTIL_ERROR_H_