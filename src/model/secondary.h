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
  Secondary() = default;
  ~Secondary() = default;
  explicit Secondary(std::initializer_list<Index> init) : data_(init) {
    VerifyRnaSize(data_.size());
  }
  explicit Secondary(std::size_t size) : data_(size, INVALID_INDEX) { VerifyRnaSize(data_.size()); }

  Secondary(Secondary&&) = default;
  Secondary& operator=(Secondary&&) = default;

  // Allow copies explicitly using the constructor.
  explicit Secondary(const Secondary&) = default;
  Secondary& operator=(const Secondary&) = delete;

  auto operator<=>(const Secondary&) const = default;

  Index& operator[](std::size_t pos) { return data_[pos]; }
  const Index& operator[](std::size_t pos) const { return data_[pos]; }

  [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
  [[nodiscard]] auto end() const noexcept { return data_.end(); }

  [[nodiscard]] auto cbegin() const noexcept { return data_.cbegin(); }
  [[nodiscard]] auto cend() const noexcept { return data_.cend(); }

  [[nodiscard]] Index size() const { return static_cast<Index>(data_.size()); }

  [[nodiscard]] bool IsPaired(int pos) const {
    assert(pos >= 0);
    return data_[pos] != INVALID_INDEX;
  }
  [[nodiscard]] bool IsUnpaired(int pos) const {
    assert(pos >= 0);
    return data_[pos] == INVALID_INDEX;
  }
  [[nodiscard]] bool IsOpeningPair(int pos) const {
    assert(pos >= 0);
    const auto pair = data_[pos];
    return pair != INVALID_INDEX && pos < pair;
  }
  [[nodiscard]] bool IsClosingPair(int pos) const {
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

  [[nodiscard]] bool PreviousPaired(int st, int en) const {
    return st != 0 && en + 1 < size() && data_[st - 1] == en + 1;
  }

  void Verify(const Primary& r) const;

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

  void Apply(Secondary& s) const {
    assert(st >= 0 && st != INVALID_INDEX);
    assert(en >= 0 && en != INVALID_INDEX);
    s[st] = en;
    s[en] = st;
  }

  void Remove(Secondary& s) const {
    assert(st >= 0 && st != INVALID_INDEX);
    assert(en >= 0 && en != INVALID_INDEX);
    s[st] = INVALID_INDEX;
    s[en] = INVALID_INDEX;
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
