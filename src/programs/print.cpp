// Copyright 2022 Eliot Courtney.
#include "programs/print.h"

#include <fmt/core.h>

#include <string>

#include "backends/common/base/model_base.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna {

namespace {

// Prints probability each pair is in a stack going inwards.
// This looks at all structures of the form *((*))* where * is any structure.
void PrintInnerStackProbs(const Primary& r, const PfnTables& pfn, const md::base::ModelBase& model,
    const erg::BoltzPseudofreeCfg& pf) {
  const auto N = As<int>(pfn.p.size());
  BoltzProbs prob(N, 0);
  if (pfn.q == ZERO_B) {
    PrintBoltzProbs(prob);
    return;
  }
  for (int i = 0; i < N; ++i)
    for (int j = i + 3; j < N; ++j) {
      const auto stacking = model.stack[r[i]][r[i + 1]][r[j - 1]][r[j]].Boltz() * pf.Paired(i, j);
      prob[i][j] = pfn.p[j][i] * pfn.p[i + 1][j - 1] * stacking / pfn.q;
    }
  PrintBoltzProbs(prob);
}

// Prints probability each pair is in a helix, where here helix is defined as a
// stack containing at least 3 adjacent pairs with no single nucleotide bulges
// between them. This looks at all structures of the form *(((*)))* where * is
// any structure.
void PrintHelixProbs(const Primary& r, const PfnTables& pfn, const md::base::ModelBase& model,
    const erg::BoltzPseudofreeCfg& pf) {
  const auto N = As<int>(pfn.p.size());
  BoltzProbs prob(N, 0);
  if (pfn.q == ZERO_B) {
    PrintBoltzProbs(prob);
    return;
  }
  for (int i = 0; i < N; ++i)
    for (int j = i + 5; j < N; ++j) {
      const auto stacking1 = model.stack[r[i]][r[i + 1]][r[j - 1]][r[j]].Boltz() * pf.Paired(i, j);
      const auto stacking2 =
          model.stack[r[i + 1]][r[i + 2]][r[j - 2]][r[j - 1]].Boltz() * pf.Paired(i + 1, j - 1);
      prob[i][j] = pfn.p[j][i] * pfn.p[i + 2][j - 2] * stacking1 * stacking2 / pfn.q;
    }
  PrintBoltzProbs(prob);
}

}  // namespace

void PrintBoltzProbs(const BoltzProbs& p) {
  const auto N = As<int>(p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      fmt::print(FLTFMT " ", p[i][j]);
    }
    fmt::print("\n");
  }
}

void PrintPfn(const BoltzSums& p) {
  const auto N = As<int>(p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      fmt::print(FLTFMT " ", p[i][j]);
    }
    fmt::print("\n");
  }
}

void PrintExtraPfnTables(const Primary& r, const PfnTables& pfn, const BackendModelPtr& model,
    const erg::PseudofreeCfg& pf) {
  pf.Verify(r);
  verify(pfn.p.size() == static_cast<std::size_t>(r.size()),
      "sequence and partition-function table must be the same length");
  const auto boltz_pf = pf.Boltz();
  std::visit(
      [&](const auto& concrete_model) {
        fmt::print("\nHelix probabilities\n");
        PrintHelixProbs(r, pfn, *concrete_model, boltz_pf);
        fmt::print("\nInner stacking probabilities\n");
        PrintInnerStackProbs(r, pfn, *concrete_model, boltz_pf);
      },
      model);
}

std::string FormatPseudofreeEnergies(std::span<const Energy> energies) {
  std::string s;
  bool first = true;
  for (const auto& e : energies) {
    if (!first) s += ",";
    first = false;
    s += fmt::format("{}", e);
  }
  return s;
}

}  // namespace mrna
