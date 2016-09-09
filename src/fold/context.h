#ifndef MEMERNA_FOLD_CONTEXT_H
#define MEMERNA_FOLD_CONTEXT_H

#include <stack>
#include "common.h"
#include "argparse.h"
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
    ZERO,
    BRUTE
  };

  static constexpr TableAlg TABLE_ALGS[] = {TableAlg::ZERO, TableAlg::ONE, TableAlg::TWO, TableAlg::THREE};
  static constexpr SuboptimalAlg SUBOPTIMAL_ALGS[] = {SuboptimalAlg::ZERO};

  context_options_t(TableAlg table_alg_ = TableAlg::ZERO,
      SuboptimalAlg suboptimal_alg_ = SuboptimalAlg::ZERO,
      energy_t subopt_delta_ = -1, int subopt_num_ = -1)
      : table_alg(table_alg_), suboptimal_alg(suboptimal_alg_),
        subopt_delta(subopt_delta_), subopt_num(subopt_num_) {}

  TableAlg table_alg;
  SuboptimalAlg suboptimal_alg;
  energy_t subopt_delta;
  int subopt_num;
};

class Context {
public:
  Context(const primary_t& r_, const energy::EnergyModelPtr em_)
      : r(r_), em(em_), options() {};
  Context(const primary_t& r_, const energy::EnergyModelPtr em_,
      context_options_t options_) : r(r_), em(em_), options(options_) {};
  Context(const Context& o) = default;

  Context() = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&& o) = delete;
  Context& operator=(Context&&) = delete;

  computed_t Fold();
  std::vector<computed_t> Suboptimal();
  const context_options_t& GetOptions() const {return options;}

private:
  const primary_t r;
  const energy::EnergyModelPtr em;
  const context_options_t options;

  void ComputeTables();
};

const std::map<std::string, ArgParse::option_t> FOLD_OPTIONS = {
    {"alg", ArgParse::option_t("which algorithm for memerna").Arg("0", {"0", "1", "2", "3", "brute"})},
    {"subopt-delta", ArgParse::option_t("maximum energy delta from minimum").Arg("-1")},
    {"subopt-num", ArgParse::option_t("maximum number of reported structures").Arg("-1")}
};

context_options_t ContextOptionsFromArgParse(const ArgParse& argparse);

}
}

#endif  // MEMERNA_FOLD_CONTEXT_H
