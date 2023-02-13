// Copyright 2021 E.
#ifndef UTIL_ERROR_H_
#define UTIL_ERROR_H_

#include <fmt/core.h>

#include <stdexcept>

#ifdef USE_BOOST
#include <boost/stacktrace.hpp>
// Like assert, but can't be disabled.
#define verify(expr, ...)                                                   \
  do {                                                                      \
    if (!(expr)) [[unlikely]] {                                             \
      auto msg = ::fmt::format("{}:{}: ", __func__, __LINE__);              \
      msg += ::fmt::format(__VA_ARGS__);                                    \
      msg += '\n';                                                          \
      msg += boost::stacktrace::to_string(boost::stacktrace::stacktrace()); \
      throw ::std::runtime_error(msg);                                      \
    }                                                                       \
  } while (0)
#else
// Like assert, but can't be disabled.
#define verify(expr, ...)                                      \
  do {                                                         \
    if (!(expr)) [[unlikely]] {                                \
      auto msg = ::fmt::format("{}:{}: ", __func__, __LINE__); \
      msg += ::fmt::format(__VA_ARGS__);                       \
      msg += '\n';                                             \
      throw ::std::runtime_error(msg);                         \
    }                                                          \
  } while (0)
#endif

#define fatal(...)              \
  do {                          \
    verify(false, __VA_ARGS__); \
    __builtin_unreachable();    \
  } while (0)

#define bug() fatal("bug")

#endif  // UTIL_ERROR_H_
