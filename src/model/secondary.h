// Copyright 2022 Eliot Courtney.
#ifndef MODEL_SECONDARY_H_
#define MODEL_SECONDARY_H_

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <string>
#include <tuple>
#include <vector>

#include "model/primary.h"

namespace mrna {

// Stores a secondary structure as a vector of indices. The index at a position
// is the position of the base it is paired with. If the base is not paired, the
// index is -1.
class Secondary {
 public:
  constexpr Secondary() = default;
  constexpr ~Secondary() = default;
  constexpr explicit Secondary(std::initializer_list<int> init) : data_(init) {}
  constexpr explicit Secondary(std::size_t size) : data_(size, -1) {}

  constexpr Secondary(Secondary&&) = default;
  constexpr Secondary& operator=(Secondary&&) = default;

  // Allow copies explicitly using the constructor.
  constexpr explicit Secondary(const Secondary&) = default;
  Secondary& operator=(const Secondary&) = delete;

  constexpr auto operator<=>(const Secondary&) const = default;

  constexpr int& operator[](std::size_t pos) { return data_[pos]; }
  constexpr const int& operator[](std::size_t pos) const { return data_[pos]; }

  [[nodiscard]] constexpr auto begin() const noexcept { return data_.begin(); }
  [[nodiscard]] constexpr auto end() const noexcept { return data_.end(); }

  [[nodiscard]] constexpr auto cbegin() const noexcept { return data_.cbegin(); }
  [[nodiscard]] constexpr auto cend() const noexcept { return data_.cend(); }

  [[nodiscard]] constexpr std::size_t size() const { return data_.size(); }

  void reset(std::size_t size) {
    data_.resize(size);
    std::fill(data_.begin(), data_.end(), -1);
  }

  [[nodiscard]] inline constexpr bool PreviousPaired(int st, int en) const {
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

  constexpr void Apply(Secondary& s) const {
    assert(st >= 0);
    assert(en >= 0);
    s[st] = en;
    s[en] = st;
  }

  constexpr void MaybeApply(Secondary& s) const {
    if (IsValid()) Apply(s);
  }
};

std::tuple<Primary, Secondary> ParseSeqDb(
    const std::string& prim_str, const std::string& pairs_str);

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
