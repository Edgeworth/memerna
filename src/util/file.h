// Copyright 2026 Eliot Courtney.
#ifndef UTIL_FILE_H_
#define UTIL_FILE_H_

#include <fstream>
#include <string>

#include "util/error.h"

namespace mrna {

[[nodiscard]] inline std::ifstream OpenFile(const std::string& path) {
  std::ifstream input(path);
  verify(input, "could not open file: {}", path);
  return input;
}

}  // namespace mrna

#endif  // UTIL_FILE_H_
