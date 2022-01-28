// Copyright 2022 Eliot Courtney.
#ifndef MODEL_CTD_H_
#define MODEL_CTD_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <string>
#include <tuple>
#include <vector>

#include "model/secondary.h"

namespace mrna {

enum Ctd : int8_t {
  // Used if CTDs are not relevant (e.g. unspecified), or if the corresponding base cannot be in a
  // CTD.
  CTD_NA,
  // Used if the corresponding base is explicitly not in a CTD, but could potentially be.
  CTD_UNUSED,
  CTD_3_DANGLE,
  CTD_5_DANGLE,
  CTD_MISMATCH,
  CTD_LCOAX_WITH_NEXT,
  CTD_LCOAX_WITH_PREV,
  CTD_RCOAX_WITH_NEXT,
  CTD_RCOAX_WITH_PREV,
  CTD_FCOAX_WITH_NEXT,
  CTD_FCOAX_WITH_PREV,
  CTD_SIZE
};

class Ctds {
 public:
  Ctds() = default;
  Ctds(std::initializer_list<Ctd> init) : data_(init) {}
  explicit Ctds(std::size_t size) : data_(size, CTD_NA) {}

  Ctds(Ctds&&) = default;
  Ctds& operator=(Ctds&&) = default;

  // Allow copies explicitly using the constructor.
  explicit Ctds(const Ctds&) = default;
  Ctds& operator=(const Ctds&) = delete;

  auto operator<=>(const Ctds&) const = default;

  Ctd& operator[](std::size_t pos) { return data_[pos]; }
  const Ctd& operator[](std::size_t pos) const { return data_[pos]; }

  auto begin() const noexcept { return data_.begin(); }
  auto end() const noexcept { return data_.end(); }

  auto cbegin() const noexcept { return data_.cbegin(); }
  auto cend() const noexcept { return data_.cend(); }

  std::size_t size() const { return data_.size(); }

  void reset(std::size_t size) {
    data_.resize(size);
    std::fill(data_.begin(), data_.end(), CTD_NA);
  }

  std::string ToString(const Secondary& s) const;
  static bool IsCtdString(const std::string& ctd_str);

 private:
  std::vector<Ctd> data_;
};

std::tuple<Primary, Secondary, Ctds> ParsePrimaryCtdString(
    const std::string& prim_str, const std::string& ctd_str);

}  // namespace mrna

#endif  // MODEL_CTD_H_
