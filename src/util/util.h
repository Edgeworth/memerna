#ifndef UTIL_UTIL_H_
#define UTIL_UTIL_H_

namespace mrna {

template <class... T>
struct overloaded : T... {
  using T::operator()...;
};

template <class... T>
overloaded(T...) -> overloaded<T...>;

}  // namespace mrna

#endif  // UTIL_UTIL_H_
