#ifndef MEMERNA_MINMAXHEAP_H
#define MEMERNA_MINMAXHEAP_H

#include <algorithm>
#include "common.h"

namespace memerna {

template<typename T, typename Comparator = std::less<T>>
class MinMaxHeap {
public:
  MinMaxHeap(const Comparator& cmp_ = Comparator()) : data(1), cmp(cmp_) {}

  const T& Min() const {
    assert(data.size() > 1u);
    if (data.size() > 3u)
      return std::min(data[2], data[3], cmp);
    return data[data.size() - 1];
  }

  const T& Max() const {
    assert(data.size() > 1u);
    return data[1];
  }

  std::size_t Size() const {
    return data.size() - 1;
  }

  void Insert(const T& val) {
    unsigned int idx = static_cast<unsigned int>(data.size());
    bool is_max = __builtin_clz(idx) % 2;
    data.push_back(val);
    if (is_max) {
      if (idx != 1 && cmp(data[idx], data[idx / 2])) {
        std::swap(data[idx], data[idx / 2]);
        SiftUpMin(idx / 2);
      } else {
        SiftUpMax(idx);
      }
    } else {
      if (idx != 1 && cmp(data[idx / 2], data[idx])) {
        std::swap(data[idx], data[idx / 2]);
        SiftUpMax(idx / 2);
      } else {
        SiftUpMin(idx);
      }
    }
  }

  void PopMin() {
    assert(data.size() > 1u);
    int min_idx = 1;
    if (data.size() > 2u && cmp(data[2], data[1])) min_idx = 2;
    if (data.size() > 3u && cmp(data[3], data[2])) min_idx = 3;
    std::swap(data[min_idx], data.back());
    data.pop_back();
    // Swapped in element should always be <= the root, since it already existed.
    assert(!cmp(data[1], data[min_idx]));
    SiftDownMin(min_idx);
  }

  void PopMax() {
    assert(data.size() > 1u);
    std::swap(data[1], data.back());
    data.pop_back();
    SiftDownMax(1);
  }

  T* begin() {return &data[1];}
  T* end() {return &data[data.size()];}

private:
  void SiftUpMin(int idx) {
    while (idx / 4 != 0) {
      if (cmp(data[idx], data[idx / 4])) std::swap(data[idx], data[idx / 4]);
      else return;
      idx /= 4;
    }
  }

  void SiftUpMax(int idx) {
    while (idx / 4 != 0) {
      if (cmp(data[idx / 4], data[idx])) std::swap(data[idx], data[idx / 4]);
      else return;
      idx /= 4;
    }
  }

  void SiftDownMax(int idx) {
    const int N = int(data.size());
    // While we have at least one child.
    while (2 * idx < N) {
      // Maintain direct child constraints.
      if (cmp(data[idx], data[2 * idx])) std::swap(data[idx], data[2 * idx]);
      if (2 * idx + 1 < N && cmp(data[idx], data[2 * idx + 1])) std::swap(data[idx], data[2 * idx + 1]);
      // Maintain grandchild constraints.
      if (4 * idx >= N) return;
      int max_idx = 4 * idx;
      if (4 * idx + 1 < N && cmp(data[max_idx], data[4 * idx + 1])) max_idx = 4 * idx + 1;
      if (4 * idx + 2 < N && cmp(data[max_idx], data[4 * idx + 2])) max_idx = 4 * idx + 2;
      if (4 * idx + 3 < N && cmp(data[max_idx], data[4 * idx + 3])) max_idx = 4 * idx + 3;
      // If our max grandchild is not bigger than us, we are done.
      if (!cmp(data[idx], data[max_idx])) return;
      std::swap(data[idx], data[max_idx]);
      idx = max_idx;
    }
  }

  void SiftDownMin(int idx) {
    const int N = int(data.size());
    // While we have at least one child.
    while (2 * idx < N) {
      // Maintain direct child constraints.
      if (cmp(data[2 * idx], data[idx])) std::swap(data[idx], data[2 * idx]);
      if (2 * idx + 1 < N && cmp(data[2 * idx + 1], data[idx])) std::swap(data[idx], data[2 * idx + 1]);
      // Maintain grandchild constraints.
      if (4 * idx >= N) return;
      int min_idx = 4 * idx;
      if (4 * idx + 1 < N && cmp(data[4 * idx + 1], data[min_idx])) min_idx = 4 * idx + 1;
      if (4 * idx + 2 < N && cmp(data[4 * idx + 2], data[min_idx])) min_idx = 4 * idx + 2;
      if (4 * idx + 3 < N && cmp(data[4 * idx + 3], data[min_idx])) min_idx = 4 * idx + 3;
      // If our min grandchild is not smaller than us, then stop.
      if (!cmp(data[min_idx], data[idx])) return;
      std::swap(data[idx], data[min_idx]);
      idx = min_idx;
    }
  }

  std::vector<T> data;
  Comparator cmp;
};

}

#endif  // MEMERNA_MINMAXHEAP_H
