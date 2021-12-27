// Copyright 2021 Eliot Courtney.
#ifndef UTIL_MACROS_H_
#define UTIL_MACROS_H_

#include <cstdio>

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

#define unimplemented() verify(false, "unimplemented")

#define bug() verify(false, "bug")

#endif  // UTIL_MACROS_H_
