#ifndef MEMERNA_FOLD_H
#define MEMERNA_FOLD_H

#include <stack>
#include "argparse.h"
#include "array.h"
#include "common.h"
#include "energy/energy.h"
#include "fold/fold_internal.h"

namespace memerna {
namespace fold {

// DP arrays
enum {
  DP_P,  // For the paired array.
  DP_U,  // For the unpaired array.
  DP_U2, // Contains at least two branches.
  DP_U_WC,  // Unpaired but must start with a branch not involved in a CTD interaction that is not GU.
  DP_U_GU,  // Unpaired but must start with a branch not involved in a CTD interaction that is GU.
  DP_U_RCOAX,  // Unpaired but must start with a branch involved in a right coaxial stack - includes energy for it.
  DP_SIZE
};

enum {
  EXT,
  EXT_WC,  // Must start with a branch not involved in an interaction that is Watson-Crick
  EXT_GU,  // Must start with a branch not involved in an interaction that is GU
  EXT_RCOAX,  // Must start with a branch, that branch is involved backwards in a right coaxial stack.
  EXT_SIZE
};

struct index_t {
  int st, en, a;

  index_t() = default;
  index_t(int st_, int en_, int a_) : st(st_), en(en_), a(a_) {}

  bool operator==(const index_t& o) const {
    return st == o.st && en == o.en && a == o.a;
  }

  bool operator<(const index_t& o) const {
    if (st != o.st) return st < o.st;
    if (en != o.en) return en < o.en;
    if (a != o.a) return a < o.a;
    return false;
  }
};

struct context_options_t {
  enum class TableAlg {
    ZERO,
    ONE,
    TWO,
    THREE
  };

  enum class SuboptimalAlg {
    ZERO
  };

  static constexpr TableAlg TABLE_ALGS[] = {TableAlg::ZERO, TableAlg::ONE, TableAlg::TWO, TableAlg::THREE};
  static constexpr SuboptimalAlg SUBOPTIMAL_ALGS[] = {SuboptimalAlg::ZERO};

  context_options_t(TableAlg table_alg_ = TableAlg::ZERO,
      SuboptimalAlg suboptimal_alg_ = SuboptimalAlg::ZERO,
      energy_t subopt_energy_ = -1, int subopt_num_ = -1)
      : table_alg(table_alg_), suboptimal_alg(suboptimal_alg_),
        subopt_energy(subopt_energy_), subopt_num(subopt_num_) {}

  TableAlg table_alg;
  SuboptimalAlg suboptimal_alg;
  energy_t subopt_energy;
  int subopt_num;
};

class Context {
public:
  Context(const primary_t& r_, const energy::EnergyModel& em_);
  Context(const primary_t& r_, const energy::EnergyModel& em_, context_options_t options_);

  Context() = delete;
  Context(const Context&) = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&& o) : r(o.r), em(o.em), options(o.options), N(o.N), pc(o.pc), tables_computed(o.tables_computed) {
    arr = std::move(o.arr);
    exterior = std::move(o.exterior);
    o.tables_computed = false;
  }
  Context& operator=(Context&&) = delete;

  computed_t Fold();
  std::vector<computed_t> Suboptimal();

  const array3d_t<energy_t, DP_SIZE>& GetDpState() const {return arr;};
  const array2d_t<energy_t, EXT_SIZE>& GetExteriorState() const {return exterior;};
  const energy::EnergyModel& GetEnergyModel() const {return em;};
  const primary_t& GetPrimary() const {return r;}
  const context_options_t& GetOptions() const {return options;}

  energy_t FastTwoLoop(int ost, int oen, int ist, int ien);
  energy_t FastHairpin(int st, int en);

private:
  const primary_t r;
  const energy::EnergyModel em;
  const context_options_t options;
  const int N;
  const internal::precomp_t pc;

  bool tables_computed;
  array3d_t<energy_t, DP_SIZE> arr;
  array2d_t<energy_t, EXT_SIZE> exterior;

  void ComputeTables0();
  void ComputeTables1();
  void ComputeTables2();
  void ComputeTables3();
  void ComputeExterior();
  void EnsureComputed();

  computed_t Traceback();
};

const std::map<std::string, ArgParse::option_t> FOLD_OPTIONS = {
    {"alg", ArgParse::option_t("which algorithm for memerna").Arg("0", {"0", "1", "2", "3"})},
    {"subopt-delta", ArgParse::option_t("maximum energy delta from minimum").Arg("-1")},
    {"subopt-num", ArgParse::option_t("maximum number of reported structures").Arg("-1")}
};

context_options_t ContextOptionsFromArgParse(const ArgParse& argparse);

}
}

#endif //MEMERNA_FOLD_H
