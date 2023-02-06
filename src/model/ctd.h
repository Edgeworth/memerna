// Copyright 2022 Eliot Courtney.
#ifndef MODEL_CTD_H_
#define MODEL_CTD_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna {

enum Ctd : uint8_t {
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
  CTD_RC_WITH_NEXT,
  CTD_RC_WITH_PREV,
  CTD_FCOAX_WITH_NEXT,
  CTD_FCOAX_WITH_PREV,
  CTD_SIZE
};

const char* CtdToName(Ctd ctd);

// Branch representation of CTDs:
//
// A CTDs used with reference to a list of branch indices. Depending on which
// side of the branch's index is used, it refers to the branch as an outer loop
// or an inner loop. By convention, if there is an outer loop, it comes first
// (not last), but it does not matter to these functions. Outer loops are
// represented by branch indices referring to the right side of the loop (i.e.
// s[branch] < branch).
using BranchCtd = std::deque<std::pair<Ctd, Energy>>;

// Per-base representation of CTDs:
//
// An array the same length as the sequence, which contains CTD identifiers.
// For each base-pair at a branch, we store a CTD identifier or nothing.
// Since there are two locations per base-pair, we can store up to two things.
// At the index of the right side of a branch, we store the CTD identifier
// for this branch for when it is an outer loop (i.e. closing a multiloop).
// At the index of the left side, we store the CTD identifier for this branch
// for when it is inside a multiloop. This is because in this situation
// ((...).(...))(...), the first loop can coaxially stack twice, one when it is
// an outer loop, and the other with the loop on the far right.
class Ctds {
 public:
  Ctds() = default;
  ~Ctds() = default;
  explicit Ctds(std::initializer_list<Ctd> init) : data_(init) {}
  explicit Ctds(std::size_t size) : data_(size, CTD_NA) {}

  Ctds(Ctds&&) = default;
  Ctds& operator=(Ctds&&) = default;

  // Allow copies explicitly using the constructor.
  explicit Ctds(const Ctds&) = default;
  Ctds& operator=(const Ctds&) = delete;

  constexpr auto operator<=>(const Ctds&) const = default;

  Ctd& operator[](std::size_t pos) { return data_[pos]; }
  const Ctd& operator[](std::size_t pos) const { return data_[pos]; }

  [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
  [[nodiscard]] auto end() const noexcept { return data_.end(); }

  [[nodiscard]] auto cbegin() const noexcept { return data_.cbegin(); }
  [[nodiscard]] auto cend() const noexcept { return data_.cend(); }

  [[nodiscard]] std::size_t size() const { return data_.size(); }

  void reset(std::size_t size) {
    data_.resize(size);
    std::fill(data_.begin(), data_.end(), CTD_NA);
  }

  [[nodiscard]] std::string ToString(const Secondary& s) const;
  static bool IsCtdString(const std::string& ctd_str);

 private:
  std::vector<Ctd> data_;
};

std::tuple<Primary, Secondary, Ctds> ParseSeqCtdString(
    const std::string& prim_str, const std::string& ctd_str);

}  // namespace mrna

#endif  // MODEL_CTD_H_
