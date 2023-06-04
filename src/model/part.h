// Copyright 2023 Eliot Courtney.
#ifndef MODEL_PART_H_
#define MODEL_PART_H_

#include <utility>

#include "model/energy.h"
#include "util/array.h"

namespace mrna {

// Holds the Boltzmann sums used in the partition function.
using BoltzSums = Array2D<BoltzEnergy>;

using BoltzProbs = Array2D<BoltzEnergy>;

struct Part {
 public:
  Part() = default;
  Part(BoltzSums p, BoltzEnergy q, BoltzProbs prob)
      : p(std::move(p)), q(q), prob(std::move(prob)) {}
  Part(BoltzSums p_, BoltzEnergy q_) : p(std::move(p_)), q(q_), prob(Prob(p, q)) {}

  [[nodiscard]] static BoltzProbs Prob(const BoltzSums& p, const BoltzEnergy& q) {
    const int N = static_cast<int>(p.size());
    BoltzProbs prob(N, 0);
    for (int i = 0; i < N; ++i)
      for (int j = i; j < N; ++j) prob[i][j] = p[i][j] * p[j][i] / q;
    return prob;
  }

  void RecomputeProb() { prob = Prob(p, q); }

  BoltzSums p;
  BoltzEnergy q{};
  BoltzProbs prob{};
};

}  // namespace mrna

#endif  // MODEL_PART_H_
