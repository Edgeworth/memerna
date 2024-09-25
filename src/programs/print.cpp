// Copyright 2022 Eliot Courtney.
#include "programs/print.h"

#include <fmt/core.h>

#include "backends/common/base/model_base.h"

namespace mrna {

void PrintBoltzProbs(const BoltzProbs& p) {
  const int N = static_cast<int>(p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      fmt::print(FLTFMT " ", p[i][j]);
    }
    fmt::print("\n");
  }
}

void PrintPfn(const BoltzSums& p) {
  const int N = static_cast<int>(p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      fmt::print(FLTFMT " ", p[i][j]);
    }
    fmt::print("\n");
  }
}

// Prints probability each pair is in a stack going inwards.
// This looks at all structures of the form *((*))* where * is any structure.
void PrintInnerStackProbs(
    const Primary& r, const PfnTables& pfn, const md::base::ModelBase& model) {
  const int N = static_cast<int>(pfn.p.size());
  BoltzProbs prob(N, 0);
  for (int i = 0; i < N; ++i)
    for (int j = i + 3; j < N; ++j) {
      const auto stacking = model.stack[r[i]][r[i + 1]][r[j - 1]][r[j]].Boltz();
      prob[i][j] = pfn.p[j][i] * pfn.p[i + 1][j - 1] * stacking / pfn.q;
    }
  PrintBoltzProbs(prob);
}

// Prints probability each pair is in a helix, where here helix is defined as a stack containing at
// least 3 adjacent pairs with no single nucleotide bulges between them. This looks at all
// structures of the form *(((*)))* where * is any structure.
void PrintHelixProbs(const Primary& r, const PfnTables& pfn, const md::base::ModelBase& model) {
  const int N = static_cast<int>(pfn.p.size());
  BoltzProbs prob(N, 0);
  for (int i = 0; i < N; ++i)
    for (int j = i + 5; j < N; ++j) {
      const auto stacking1 = model.stack[r[i]][r[i + 1]][r[j - 1]][r[j]].Boltz();
      const auto stacking2 = model.stack[r[i + 1]][r[i + 2]][r[j - 2]][r[j - 1]].Boltz();
      prob[i][j] = pfn.p[j][i] * pfn.p[i + 2][j - 2] * stacking1 * stacking2 / pfn.q;
    }
  PrintBoltzProbs(prob);
}

}  // namespace mrna
