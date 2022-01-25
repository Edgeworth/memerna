// Copyright 2022 E.
#ifndef MODEL_PRIMARY_H_
#define MODEL_PRIMARY_H_

#include <cstddef>
#include <string>
#include <vector>

#include "model/base.h"

namespace mrna {

class Primary {
 public:
  Primary() = default;
  explicit Primary(std::size_t size) : data_(size, A) {}

  Primary(Primary&&) = default;
  Primary& operator=(Primary&&) = default;

  // Allow copies explicitly using the constructor.
  explicit Primary(const Primary&) = default;
  Primary& operator=(const Primary&) = delete;

  auto operator<=>(const Primary&) const = default;

  Base& operator[](std::size_t pos) { return data_[pos]; }
  const Base& operator[](std::size_t pos) const { return data_[pos]; }

  auto begin() const noexcept { return data_.begin(); }
  auto end() const noexcept { return data_.end(); }

  auto cbegin() const noexcept { return data_.cbegin(); }
  auto cend() const noexcept { return data_.cend(); }

  std::size_t size() const { return data_.size(); }

  std::string ToString() const;

  static Primary Random(int length);
  static Primary FromString(const std::string& s);

 private:
  std::vector<Base> data_;
};

}  // namespace mrna

#endif  // MODEL_PRIMARY_H_
