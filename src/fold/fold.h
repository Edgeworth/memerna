#ifndef MEMERNA_FOLD_H
#define MEMERNA_FOLD_H

#include "argparse.h"
#include "array.h"
#include "constants.h"
#include "common.h"
#include "energy/energy.h"

namespace memerna {
namespace fold {

const int MAX_SPECIAL_HAIRPIN_SZ = 6;

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

// Split candidates up into several lists.
// In general, for each array we need a new candidate list (except for U and U2 which mirror each other very
// closely). We also need another candidate list for forward RCOAX since we can't use its energy value directly, not
// knowing it. Same with flush coaxial stacks.
enum {
  CAND_P_MISMATCH,  // No monotonicity.
  CAND_P_OUTER,  // No monotonicity.
  CAND_P_FLUSH,  // No monotonicity.
  CAND_U,
  CAND_U_LCOAX,  // No monotonicity.
  CAND_U_RCOAX_FWD,  // No monotonicity.
  CAND_U_WC_FLUSH,  // No monotonicity.
  CAND_U_GU_FLUSH,  // No monotonicity.
  CAND_U_WC,
  CAND_U_GU,
  CAND_U_RCOAX,
  CAND_SIZE
};

enum {
  CAND_EN_P_MISMATCH,
  CAND_EN_P_OUTER,
  CAND_EN_P_FLUSH,
  CAND_EN_SIZE
};

struct cand_t {
  energy_t energy;
  int idx;
};

struct hairpin_precomp_t {
  hairpin_precomp_t() : num_c(0) {
    memset(special, constants::MAX_E & 0xFF, sizeof(special));
  }

  energy_t special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};


int MaxNumContiguous(const primary_t& primary);
energy_t FastTwoLoop(int ost, int oen, int ist, int ien);
std::vector<hairpin_precomp_t> PrecomputeFastHairpin();
energy_t FastHairpin(int st, int en, const std::vector<hairpin_precomp_t>& precomp);

inline bool IsNotLonely(int st, int en) {
  return (en - st - 3 >= constants::HAIRPIN_MIN_SZ && CanPair(r[st + 1], r[en - 1])) ||
      (st > 0 && en < int(r.size() - 1) && CanPair(r[st - 1], r[en + 1]));
}


struct fold_state_t {
  array3d_t<energy_t, DP_SIZE> dp_table;
  array2d_t<energy_t, EXT_SIZE> ext_table;
};

array3d_t<energy_t, DP_SIZE> ComputeTables3();
array3d_t<energy_t, DP_SIZE> ComputeTables2();
array3d_t<energy_t, DP_SIZE> ComputeTables1();
array3d_t<energy_t, DP_SIZE> ComputeTables0();

void InitFold();
computed_t Fold3(const primary_t& primary, fold_state_t* fold_state = nullptr);
computed_t Fold2(const primary_t& primary, fold_state_t* fold_state = nullptr);
computed_t Fold1(const primary_t& primary, fold_state_t* fold_state = nullptr);
computed_t Fold0(const primary_t& primary, fold_state_t* fold_state = nullptr);
computed_t FoldBruteForce(const primary_t& primary, fold_state_t* fold_state = nullptr);

typedef computed_t (fold_fn_t)(const primary_t&, fold_state_t* fold_state);

fold_fn_t* const FOLD_FUNCTIONS[] = {&Fold0, &Fold1, &Fold2, &Fold3};

const std::map<std::string, ArgParse::option_t> FOLD_OPTIONS = {
    {"alg", ArgParse::option_t("which algorithm for memerna").Arg("0", {"0", "1", "2", "3", "brute"})}
};

fold::fold_fn_t* FoldFunctionFromArgParse(const ArgParse& argparse);

}
}

#endif //MEMERNA_FOLD_H
