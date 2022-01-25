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

class Secondary {
 public:
  Secondary() = default;
  Secondary(std::initializer_list<int> init) : data_(init) {}
  explicit Secondary(std::size_t size) : data_(size, -1) {}

  Secondary(Secondary&&) = default;
  Secondary& operator=(Secondary&&) = default;

  // Allow copies explicitly using the constructor.
  explicit Secondary(const Secondary&) = default;
  Secondary& operator=(const Secondary&) = delete;

  auto operator<=>(const Secondary&) const = default;

  int& operator[](std::size_t pos) { return data_[pos]; }
  const int& operator[](std::size_t pos) const { return data_[pos]; }

  auto begin() const noexcept { return data_.begin(); }
  auto end() const noexcept { return data_.end(); }

  auto cbegin() const noexcept { return data_.cbegin(); }
  auto cend() const noexcept { return data_.cend(); }

  std::size_t size() const { return data_.size(); }

  void reset(std::size_t size) {
    data_.resize(size);
    std::fill(data_.begin(), data_.end(), -1);
  }

  static Secondary FromDotBracket(const std::string& pairs_str);
  std::string ToDotBracket() const;

 private:
  std::vector<int> data_;
};

std::tuple<Primary, Secondary> ParsePrimaryDotBracket(
    const std::string& prim_str, const std::string& pairs_str);

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
