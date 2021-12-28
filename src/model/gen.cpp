// Copyright 2021 E.
#include "model/gen.h"

namespace mrna {

Primary GenerateRandomPrimary(int length) {
  static thread_local std::mt19937 eng;
  std::uniform_int_distribution<int> dist(0, 3);
  Primary r(std::size_t(length), 0);
  for (int i = 0; i < length; ++i) r[i] = Base(dist(eng));
  return r;
}

}  // namespace mrna
