// Copyright 2022 Eliot Courtney.
#ifndef MODEL_PRIMARY_H_
#define MODEL_PRIMARY_H_

#include <cstddef>
#include <string>
#include <vector>

#include "model/base.h"

namespace mrna {

class Primary {
 public:
  constexpr Primary() = default;
  constexpr ~Primary() = default;
  constexpr explicit Primary(std::size_t size) : data_(size, A) {}

  constexpr Primary(Primary&&) = default;
  constexpr Primary& operator=(Primary&&) = default;

  // Allow copies explicitly using the constructor.
  constexpr explicit Primary(const Primary&) = default;
  Primary& operator=(const Primary&) = delete;

  constexpr auto operator<=>(const Primary&) const = default;

  constexpr Base& operator[](std::size_t pos) { return data_[pos]; }
  constexpr const Base& operator[](std::size_t pos) const { return data_[pos]; }

  [[nodiscard]] constexpr auto begin() const noexcept { return data_.begin(); }
  [[nodiscard]] constexpr auto end() const noexcept { return data_.end(); }

  [[nodiscard]] constexpr auto cbegin() const noexcept { return data_.cbegin(); }
  [[nodiscard]] constexpr auto cend() const noexcept { return data_.cend(); }

  [[nodiscard]] constexpr std::size_t size() const { return data_.size(); }

  [[nodiscard]] std::string ToSeq() const;

  // Treat the sequence as a number and perform an increment.
  void Increment();

  static Primary Random(int length);
  static Primary FromSeq(const std::string& s);

 private:
  std::vector<Base> data_;
};

}  // namespace mrna

#endif  // MODEL_PRIMARY_H_
