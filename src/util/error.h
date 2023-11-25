// Copyright 2021 Eliot Courtney.
#ifndef UTIL_ERROR_H_
#define UTIL_ERROR_H_

#include <fmt/format.h>

#include <stdexcept>

#define verify(expr, ...)                                      \
  do {                                                         \
    if (!(expr)) [[unlikely]] {                                \
      auto msg = ::fmt::format("{}:{}: ", __func__, __LINE__); \
      msg += ::fmt::format(__VA_ARGS__);                       \
      msg += '\n';                                             \
      throw ::std::runtime_error(msg);                         \
    }                                                          \
  } while (0) /* NOLINT(cppcoreguidelines-avoid-do-while) */

#define fatal(...)              \
  do {                          \
    verify(false, __VA_ARGS__); \
    __builtin_unreachable();    \
  } while (0) /* NOLINT(cppcoreguidelines-avoid-do-while) */

#define bug() fatal("bug")

namespace mrna {

void InitProgram();

}  // namespace mrna

#endif  // UTIL_ERROR_H_
