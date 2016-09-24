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

std::vector<computed_t> Context::SuboptimalSorted(energy_t subopt_delta, int subopt_num) {
  std::vector<computed_t> computeds;
  Suboptimal([&computeds](const computed_t& c) { computeds.push_back(c); },
      subopt_delta, subopt_num);
  if (!IsSuboptimalSorted(subopt_delta, subopt_num))
    std::sort(computeds.begin(), computeds.end(), computed_energy_comparator_t());
  return computeds;
}

void Context::Suboptimal(std::function<void(const computed_t&)> fn,
    energy_t subopt_delta, int subopt_num) {
  if (options.suboptimal_alg == context_options_t::SuboptimalAlg::BRUTE) {
    for (const auto& computed : FoldBruteForce(r, *em, subopt_num))
      fn(computed);
    return;
  }

  ComputeTables();
  switch (options.suboptimal_alg) {
    case context_options_t::SuboptimalAlg::ZERO:
      internal::Suboptimal0(subopt_delta, subopt_num).Run(fn);
      break;
    case context_options_t::SuboptimalAlg::ONE:
      internal::Suboptimal1(subopt_delta, subopt_num).Run(fn);
      break;
    default:
      verify_expr(false, "bug - no such suboptimal algorithm %d", int(options.suboptimal_alg));
  }
}

bool Context::IsSuboptimalSorted(energy_t, int subopt_num) {
  return !(options.suboptimal_alg == context_options_t::SuboptimalAlg::ONE && subopt_num == -1);
}

}
}
