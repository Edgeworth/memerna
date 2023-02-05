// Copyright 2016 Eliot Courtney.
#ifndef API_PART_PART_H_
#define API_PART_PART_H_

#include <cassert>
#include <variant>

#include "model/energy.h"
#include "models/t04/part/part.h"

namespace mrna::part {

// Holds the Boltzmann sums used in the partition function.
using BoltzSums = Array2D<BoltzEnergy>;

using BoltzProbs = Array2D<BoltzEnergy>;

using PartState = std::variant<md::t04::PartState>;

struct Part {
  BoltzSums p;
  BoltzEnergy q;

  [[nodiscard]] BoltzProbs Prob() const {
    const int N = static_cast<int>(p.size());
    BoltzProbs prob(N, 0);
    for (int i = 0; i < N; ++i)
      for (int j = i; j < N; ++j) prob[i][j] = p[i][j] * p[j][i] / q;
    return prob;
  }
};

struct PartResult {
  PartState state;
  Part part;
  BoltzProbs prob;
};

}  // namespace mrna::part

#endif  // API_PART_PART_H_
