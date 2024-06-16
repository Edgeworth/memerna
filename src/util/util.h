// Copyright 2023 Eliot Courtney.
#ifndef UTIL_UTIL_H_
#define UTIL_UTIL_H_

#include <algorithm>
#include <cassert>
#include <type_traits>

namespace mrna {

template <class... T>
struct overloaded : T... {
  using T::operator()...;
};

template <class... T>
overloaded(T...) -> overloaded<T...>;

template <typename T>
constexpr auto Decay(T& a) {
  return reinterpret_cast<std::remove_all_extents_t<T>*>(&a);
}

// Only works with [0, 2N).
inline int FastMod(int a, int m) {
  assert(a < 2 * m);
  if (a >= m) return a - m;
  return a;
}

template <typename Container, typename Value>
inline bool Contains(const Container& c, const Value& v) {
  return std::find(c.begin(), c.end(), v) != c.end();
}

}  // namespace mrna

#endif  // UTIL_UTIL_H_
