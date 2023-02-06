// Copyright 2016 Eliot Courtney.
#ifndef API_PART_H_
#define API_PART_H_

#include <variant>

#include "model/part.h"
#include "models/t04/part/part.h"

namespace mrna::part {

using PartState = std::variant<std::monostate, md::t04::PartState>;

struct PartResult {
  PartState state;
  Part part;
};

}  // namespace mrna::part

#endif  // API_PART_H_
