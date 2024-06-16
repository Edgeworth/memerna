// Copyright 2021 Eliot Courtney.
#ifndef UTIL_ERROR_H_
#define UTIL_ERROR_H_

#include <fmt/format.h>

#define funcname() __PRETTY_FUNCTION__

#define verify(expr, ...)                                        \
  do {                                                           \
    if (!(expr)) [[unlikely]] {                                  \
      auto msg = ::fmt::format("{}:{}: ", funcname(), __LINE__); \
      msg += ::fmt::format(__VA_ARGS__);                         \
      msg += '\n';                                               \
      throw ::std::runtime_error(msg);                           \
    }                                                            \
  } while (0)

#define fatal(...)              \
  do {                          \
    verify(false, __VA_ARGS__); \
    __builtin_unreachable();    \
  } while (0)

#ifdef NDEBUG
#define unreachable() __builtin_unreachable()
#else

#define unreachable()     \
  do {                    \
    fatal("unreachable"); \
  } while (0)

#endif

namespace mrna {

void InitProgram();

}  // namespace mrna

#endif  // UTIL_ERROR_H_
