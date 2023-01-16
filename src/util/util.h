// Copyright 2023 Eliot Courtney.
#ifndef UTIL_UTIL_H_
#define UTIL_UTIL_H_

#include <type_traits>

namespace mrna {

template <class... T>
struct overloaded : T... {
  using T::operator()...;
};

template <class... T>
overloaded(T...) -> overloaded<T...>;

template <typename T>
constexpr auto Decay(T& a) {  // NOLINT
  return reinterpret_cast<std::remove_all_extents_t<T>*>(&a);
}

}  // namespace mrna

#endif  // UTIL_UTIL_H_
