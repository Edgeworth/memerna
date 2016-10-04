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
#ifndef MEMERNA_FOLD_CONTEXT_H
#define MEMERNA_FOLD_CONTEXT_H

#include <stack>
#include "argparse.h"
#include "common.h"
#include "energy/energy.h"
#include "fold/fold.h"

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
  std::vector<computed_t> SuboptimalIntoVector(bool sorted,
      energy_t subopt_delta = -1, int subopt_num = -1);
  int Suboptimal(SuboptimalCallback fn, bool sorted,
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
