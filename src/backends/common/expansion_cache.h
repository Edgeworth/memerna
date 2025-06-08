// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_COMMON_EXPANSION_CACHE_H_
#define BACKENDS_COMMON_EXPANSION_CACHE_H_

#include <utility>
#include <vector>

#include "model/primary.h"
#include "util/lru_vector_cache.h"

namespace mrna::md {

template <typename DpIndex, typename Expansion, bool UseLru>
class ExpansionCache {
 public:
  // Use O(N^2) cache if using LRU.
  explicit ExpansionCache(const Primary& r, std::size_t max_linear_index)
      : n_(r.size()), cache_([&] {
          if constexpr (UseLru) {
            return LruVectorCache<Expansion>(max_linear_index, n_ * n_);
          } else {
            return std::vector<std::vector<Expansion>>(max_linear_index);
          }
        }()) {}

  const std::vector<Expansion>& Get(std::size_t key) {
    if constexpr (UseLru) {
      return cache_.Get(key);
    } else {
      return cache_[key];
    }
  }

  const std::vector<Expansion>& Insert(std::size_t key, std::vector<Expansion>&& value) {
    if constexpr (UseLru) {
      return cache_.Insert(key, std::move(value));
    } else {
      return cache_[key] = std::move(value);
    }
  }

 private:
  std::size_t n_;
  std::conditional_t<UseLru, LruVectorCache<Expansion>, std::vector<std::vector<Expansion>>> cache_;
};

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_EXPANSION_CACHE_H_
