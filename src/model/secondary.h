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
  Secondary() = default;
  ~Secondary() = default;
  explicit Secondary(std::initializer_list<int> init) : data_(init) {}
  explicit Secondary(std::size_t size) : data_(size, -1) {}

  Secondary(Secondary&&) = default;
  Secondary& operator=(Secondary&&) = default;

  // Allow copies explicitly using the constructor.
  explicit Secondary(const Secondary&) = default;
  Secondary& operator=(const Secondary&) = delete;

  constexpr auto operator<=>(const Secondary&) const = default;

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

  static Secondary FromDb(const std::string& pairs_str);  // Dotbracket
  [[nodiscard]] std::string ToDb() const;

 private:
  std::vector<int> data_;
};

std::tuple<Primary, Secondary> ParseSeqDb(
    const std::string& prim_str, const std::string& pairs_str);

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
