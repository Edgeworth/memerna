#ifndef MEMERNA_FOLD_CONTEXT_H
#define MEMERNA_FOLD_CONTEXT_H

#include <stack>
#include "argparse.h"
#include "common.h"
#include "energy/energy.h"

namespace memerna {
namespace fold {

struct context_options_t {
  enum class TableAlg {
    ZERO,
    ONE,
    TWO,
    THREE,
    BRUTE  // Not included in the normal table algs since exponential.
  };

  enum class SuboptimalAlg {
    ZERO, ONE, BRUTE
  };

  static constexpr TableAlg TABLE_ALGS[] = {
      TableAlg::ZERO, TableAlg::ONE, TableAlg::TWO, TableAlg::THREE};
  static constexpr SuboptimalAlg SUBOPTIMAL_ALGS[] = {SuboptimalAlg::ZERO, SuboptimalAlg::ONE};

  context_options_t(
      TableAlg table_alg_ = TableAlg::ZERO, SuboptimalAlg suboptimal_alg_ = SuboptimalAlg::ZERO)
      : table_alg(table_alg_), suboptimal_alg(suboptimal_alg_) {}

  TableAlg table_alg;
  SuboptimalAlg suboptimal_alg;
};

class Context {
public:
  Context(const primary_t& r_, const energy::EnergyModelPtr em_) : r(r_), em(em_), options() {
    verify_expr(r.size() > 0u, "cannot fold zero length RNA");
  };
  Context(const primary_t& r_, const energy::EnergyModelPtr em_, context_options_t options_)
      : r(r_), em(em_), options(options_) {
    verify_expr(r.size() > 0u, "cannot fold zero length RNA");
  }
  Context(const Context& o) = default;

  Context() = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&& o) = delete;
  Context& operator=(Context&&) = delete;

  computed_t Fold();
  std::vector<computed_t> SuboptimalSorted(energy_t subopt_delta = -1, int subopt_num = -1);
  bool IsSuboptimalSorted(energy_t subopt_delta, int subopt_num);
  void Suboptimal(std::function<void(const computed_t&)> fn,
      energy_t subopt_delta = -1, int subopt_num = -1);

private:
  const primary_t r;
  const energy::EnergyModelPtr em;
  const context_options_t options;

  void ComputeTables();
};

const std::map<std::string, ArgParse::option_t> FOLD_OPTIONS = {
    {"dp-alg",
        ArgParse::option_t("which algorithm for memerna").Arg("2", {"0", "1", "2", "3", "brute"})},
    {"subopt-alg", ArgParse::option_t("which algorithm for memerna").Arg("1", {"0", "1", "brute"})}};

context_options_t ContextOptionsFromArgParse(const ArgParse& argparse);
}
}

#endif  // MEMERNA_FOLD_CONTEXT_H
