// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#include "fold/context.h"
#include "fold/brute_fold.h"
#include "fold/suboptimal0.h"
#include "fold/suboptimal1.h"

namespace memerna {
namespace fold {

using namespace energy;

constexpr context_options_t::TableAlg context_options_t::TABLE_ALGS[];
constexpr context_options_t::SuboptimalAlg context_options_t::SUBOPTIMAL_ALGS[];

context_options_t ContextOptionsFromArgParse(const ArgParse& argparse) {
  context_options_t options;
  auto dp_alg = argparse.GetOption("dp-alg");
  if (dp_alg == "0") {
    options.table_alg = context_options_t::TableAlg::ZERO;
  } else if (dp_alg == "1") {
    options.table_alg = context_options_t::TableAlg::ONE;
  } else if (dp_alg == "2") {
    options.table_alg = context_options_t::TableAlg::TWO;
  } else if (dp_alg == "3") {
    options.table_alg = context_options_t::TableAlg::THREE;
  } else if (dp_alg == "brute") {
    options.table_alg = context_options_t::TableAlg::BRUTE;
  } else {
    verify_expr(false, "unknown fold option");
  }
  auto subopt_alg = argparse.GetOption("subopt-alg");
  if (subopt_alg == "0") {
    options.suboptimal_alg = context_options_t::SuboptimalAlg::ZERO;
  } else if (subopt_alg == "1") {
    options.suboptimal_alg = context_options_t::SuboptimalAlg::ONE;
  } else if (subopt_alg == "brute") {
    options.suboptimal_alg = context_options_t::SuboptimalAlg::BRUTE;
  } else {
    verify_expr(false, "unknown fold option");
  }
  return options;
}

void Context::ComputeTables() {
  internal::SetGlobalState(r, *em);
  switch (options.table_alg) {
    case context_options_t::TableAlg::ZERO:
      internal::ComputeTables0();
      break;
    case context_options_t::TableAlg::ONE:
      internal::ComputeTables1();
      break;
    case context_options_t::TableAlg::TWO:
      internal::ComputeTables2();
      break;
    case context_options_t::TableAlg::THREE:
      internal::ComputeTables3();
      break;
    default:
      verify_expr(false, "bug");
  }
  internal::ComputeExterior();
}

computed_t Context::Fold() {
  if (options.table_alg == context_options_t::TableAlg::BRUTE) return FoldBruteForce(r, *em, 1)[0];

  ComputeTables();
  internal::Traceback();
  return {{internal::gr, internal::gp}, internal::gctd, internal::genergy};
}

std::vector<computed_t> Context::SuboptimalIntoVector(bool sorted,
    energy_t subopt_delta, int subopt_num) {
  std::vector<computed_t> computeds;
  int num_structures = Suboptimal(
      [&computeds](const computed_t& c) { computeds.push_back(c); },
      sorted,
      subopt_delta, subopt_num);
  assert(num_structures == int(computeds.size()));
  return computeds;
}

int Context::Suboptimal(SuboptimalCallback fn, bool sorted,
    energy_t subopt_delta, int subopt_num) {
  if (options.suboptimal_alg == context_options_t::SuboptimalAlg::BRUTE) {
    auto computeds = FoldBruteForce(r, *em, subopt_num);
    for (const auto& computed : computeds)
      fn(computed);
    return int(computeds.size());
  }

  ComputeTables();
  switch (options.suboptimal_alg) {
    case context_options_t::SuboptimalAlg::ZERO:
      return internal::Suboptimal0(subopt_delta, subopt_num).Run(fn);
    case context_options_t::SuboptimalAlg::ONE:
      return internal::Suboptimal1(subopt_delta, subopt_num).Run(fn, sorted);
    default:
      verify_expr(false, "bug - no such suboptimal algorithm %d", int(options.suboptimal_alg));
  }
}

}
}
