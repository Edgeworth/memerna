// Copyright 2025 Eliot Courtney.
#ifndef UTIL_LRU_VECTOR_CACHE_H_
#define UTIL_LRU_VECTOR_CACHE_H_

#include <list>
#include <utility>
#include <vector>

namespace mrna {

// Somewhat special LRU cache since it holds vectors and clears the vectors when they are evicted.
template <typename Value>
class LruVectorCache {
 public:
  using Key = std::size_t;

  LruVectorCache(std::size_t max_size, std::size_t max_total)
      : lru_(), cache_(max_size, {{}, lru_.end()}), max_total_(max_total) {}

  const std::vector<Value>& Insert(Key key, std::vector<Value>&& value) {
    // Key should not already be in the cache. This means we don't need to remove it from the LRU
    // list.
    assert(cache_[key].first.empty());
    assert(key < cache_.size());
    total_ += value.size();
    lru_.push_front(key);
    cache_[key] = {std::move(value), lru_.begin()};
    MaybeEvict();
    return cache_[key].first;
  }

  const std::vector<Value>& Get(Key key) const {
    assert(key < cache_.size());
    if (cache_[key].second != lru_.end()) lru_.splice(lru_.begin(), lru_, cache_[key].second);
    return cache_[key].first;
  }

 private:
  mutable std::list<Key> lru_;
  std::vector<std::pair<std::vector<Value>, std::list<Key>::iterator>> cache_;
  std::size_t max_total_;
  std::size_t total_ = 0;

  void MaybeEvict() {
    while (total_ > max_total_) {
      auto key = lru_.back();
      total_ -= cache_[key].first.size();
      cache_[key].first.clear();
      cache_[key].first.shrink_to_fit();
      cache_[key].second = lru_.end();
      lru_.pop_back();
    }
  }
};

}  // namespace mrna

#endif  // UTIL_LRU_VECTOR_CACHE_H_
