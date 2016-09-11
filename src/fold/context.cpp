#include <stack>
#include <climits>
#include "parsing.h"
#include "fold/context.h"
#include "fold/suboptimal0.h"
#include "fold/brute_fold.h"

namespace memerna {
namespace fold {

using namespace constants;
using namespace energy;

constexpr context_options_t::TableAlg context_options_t::TABLE_ALGS[];
constexpr context_options_t::SuboptimalAlg context_options_t::SUBOPTIMAL_ALGS[];

context_options_t ContextOptionsFromArgParse(const ArgParse& argparse) {
  context_options_t options;
  auto opt = argparse.GetOption("alg");
  if (opt == "0") {
    options.table_alg = context_options_t::TableAlg::ZERO;
  } else if (opt == "1") {
    options.table_alg = context_options_t::TableAlg::ONE;
  } else if (opt == "2") {
    options.table_alg = context_options_t::TableAlg::TWO;
  } else if (opt == "3") {
    options.table_alg = context_options_t::TableAlg::THREE;
  } else if (opt == "brute") {
    options.table_alg = context_options_t::TableAlg::BRUTE;
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
  if (options.table_alg == context_options_t::TableAlg::BRUTE)
    return FoldBruteForce(r, *em, 1)[0];

  ComputeTables();
  internal::Traceback();
  return {{internal::gr, internal::gp}, internal::gctd, internal::genergy};
}

std::vector<computed_t> Context::Suboptimal(energy_t subopt_delta, int subopt_num) {
  if (options.suboptimal_alg == context_options_t::SuboptimalAlg::BRUTE)
    return FoldBruteForce(r, *em, subopt_num);

  ComputeTables();
  energy_t max_energy = internal::gext[0][internal::EXT] + subopt_delta;
  int max_structures = subopt_num;
  if (subopt_delta < 0) max_energy = constants::CAP_E;
  if (subopt_num < 0) max_structures = INT_MAX;
  switch (options.suboptimal_alg) {
    case context_options_t::SuboptimalAlg::ZERO:
      return internal::Suboptimal0(max_energy, max_structures).Run();
    default:
      verify_expr(false, "bug");
  }
}

}
}
