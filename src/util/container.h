// Copyright 2024 Eliot Courtney.
#ifndef UTIL_CONTAINER_H_
#define UTIL_CONTAINER_H_

#include <boost/container/small_vector.hpp>

namespace mrna {

template <class T, std::size_t N>
using smallvec = boost::container::small_vector<T, N>;

}  // namespace mrna

#endif  // UTIL_CONTAINER_H_
