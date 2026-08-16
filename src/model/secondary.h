// Copyright 2022 Eliot Courtney.
#ifndef MODEL_SECONDARY_H_
#define MODEL_SECONDARY_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <string>
#include <tuple>
#include <vector>

#include "model/base.h"
#include "model/primary.h"

namespace mrna {

// Stores a secondary structure as a vector of indices. The index at a position
// is the position of the base it is paired with. If the base is not paired, the
// index is -1.
class Secondary {
 public:
  Secondary() = default;
  ~Secondary() = default;
  explicit Secondary(std::initializer_list<int> init) : data_(init) {}
  explicit Secondary(std::size_t size) : data_(size, -1) {}

  Secondary(Secondary&&) = default;
  Secondary& operator=(Secondary&&) = default;

  // Allow copies explicitly using the constructor.
  explicit Secondary(const Secondary&) = default;
  Secondary& operator=(const Secondary&) = delete;

  auto operator<=>(const Secondary&) const = default;

  int& operator[](std::size_t pos) { return data_[pos]; }
  const int& operator[](std::size_t pos) const { return data_[pos]; }

  [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
  [[nodiscard]] auto end() const noexcept { return data_.end(); }

  [[nodiscard]] auto cbegin() const noexcept { return data_.cbegin(); }
  [[nodiscard]] auto cend() const noexcept { return data_.cend(); }

  [[nodiscard]] std::size_t size() const { return data_.size(); }

  void reset(std::size_t size) {
    data_.resize(size);
    std::fill(data_.begin(), data_.end(), -1);
  }

  void reset() { std::fill(data_.begin(), data_.end(), -1); }

  [[nodiscard]] bool PreviousPaired(int st, int en) const {
    return st != 0 && en != static_cast<int>(size()) - 1 && data_[st - 1] == en + 1;
  }

  static Secondary FromDb(const std::string& pairs_str);  // Dotbracket
  [[nodiscard]] std::string ToDb() const;

 private:
  std::vector<int> data_;
};

struct Pair {
  Index st = -1;
  Index en = -1;

  constexpr Pair() = default;
  constexpr Pair(int st_, int en_) : st(Index(st_)), en(Index(en_)) {
    assert(st_ == st);
    assert(en_ == en);
  }

  [[nodiscard]] constexpr bool IsValid() const { return st != -1; }

  void Apply(Secondary& s) const {
    assert(st >= 0);
    assert(en >= 0);
    s[st] = en;
    s[en] = st;
  }

  void Remove(Secondary& s) const {
    assert(st >= 0);
    assert(en >= 0);
    s[st] = -1;
    s[en] = -1;
  }

  void MaybeApply(Secondary& s) const {
    if (IsValid()) Apply(s);
  }

  void MaybeRemove(Secondary& s) const {
    if (IsValid()) Remove(s);
  }
};

std::tuple<Primary, Secondary> ParseSeqDb(
    const std::string& prim_str, const std::string& pairs_str);

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
