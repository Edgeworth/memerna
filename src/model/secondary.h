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
#include "util/util.h"

namespace mrna {

// Stores a secondary structure as a vector of indices. The index at a position
// is the position of the base it is paired with. If the base is not paired, the
// index is INVALID_INDEX.
class Secondary {
 public:
  constexpr Secondary() = default;
  constexpr ~Secondary() = default;
  explicit Secondary(std::initializer_list<Index> init) : data_(init) {
    VerifyRnaSize(data_.size());
  }
  explicit Secondary(std::size_t size) : data_(size, INVALID_INDEX) { VerifyRnaSize(data_.size()); }

  constexpr Secondary(Secondary&&) = default;
  constexpr Secondary& operator=(Secondary&&) = default;

  // Allow copies explicitly using the constructor.
  constexpr explicit Secondary(const Secondary&) = default;
  Secondary& operator=(const Secondary&) = delete;

  constexpr auto operator<=>(const Secondary&) const = default;

  constexpr Index& operator[](std::size_t pos) { return data_[pos]; }
  constexpr const Index& operator[](std::size_t pos) const { return data_[pos]; }

  [[nodiscard]] constexpr auto begin() const noexcept { return data_.begin(); }
  [[nodiscard]] constexpr auto end() const noexcept { return data_.end(); }

  [[nodiscard]] constexpr auto cbegin() const noexcept { return data_.cbegin(); }
  [[nodiscard]] constexpr auto cend() const noexcept { return data_.cend(); }

  [[nodiscard]] constexpr Index size() const { return static_cast<Index>(data_.size()); }

  [[nodiscard]] constexpr bool IsPaired(int pos) const {
    assert(pos >= 0);
    return data_[pos] != INVALID_INDEX;
  }
  [[nodiscard]] constexpr bool IsUnpaired(int pos) const {
    assert(pos >= 0);
    return data_[pos] == INVALID_INDEX;
  }
  [[nodiscard]] constexpr bool IsOpeningPair(int pos) const {
    assert(pos >= 0);
    const auto pair = data_[pos];
    return pair != INVALID_INDEX && pos < pair;
  }
  [[nodiscard]] constexpr bool IsClosingPair(int pos) const {
    assert(pos >= 0);
    const auto pair = data_[pos];
    return pair != INVALID_INDEX && pair < pos;
  }

  void reset(std::size_t size) {
    VerifyRnaSize(size);
    data_.resize(size);
    std::fill(data_.begin(), data_.end(), INVALID_INDEX);
  }

  void reset() { std::fill(data_.begin(), data_.end(), INVALID_INDEX); }

  [[nodiscard]] constexpr bool PreviousPaired(int st, int en) const {
    return st != 0 && en + 1 < size() && data_[st - 1] == en + 1;
  }

  static Secondary FromDb(const std::string& pairs_str);  // Dotbracket
  [[nodiscard]] std::string ToDb() const;

 private:
  std::vector<Index> data_;
};

struct Pair {
  Index st = INVALID_INDEX;
  Index en = INVALID_INDEX;

  constexpr Pair() = default;
  constexpr Pair(int st_, int en_) : st(As<Index>(st_)), en(As<Index>(en_)) {}

  [[nodiscard]] constexpr bool IsValid() const { return st != INVALID_INDEX; }

  constexpr void Apply(Secondary& s) const {
    assert(st >= 0 && st != INVALID_INDEX);
    assert(en >= 0 && en != INVALID_INDEX);
    s[st] = en;
    s[en] = st;
  }

  constexpr void Remove(Secondary& s) const {
    assert(st >= 0 && st != INVALID_INDEX);
    assert(en >= 0 && en != INVALID_INDEX);
    s[st] = INVALID_INDEX;
    s[en] = INVALID_INDEX;
  }

  constexpr void MaybeApply(Secondary& s) const {
    if (IsValid()) Apply(s);
  }

  constexpr void MaybeRemove(Secondary& s) const {
    if (IsValid()) Remove(s);
  }
};

std::tuple<Primary, Secondary> ParseSeqDb(
    const std::string& prim_str, const std::string& pairs_str);

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
