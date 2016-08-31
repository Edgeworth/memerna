#ifndef MEMERNA_FOLD_H
#define MEMERNA_FOLD_H

#include "argparse.h"
#include "array.h"
#include "common.h"
#include "energy/energy.h"

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

struct fold_state_t {
  array3d_t<energy_t, DP_SIZE> dp_table;
  array2d_t<energy_t, EXT_SIZE> ext_table;
};

computed_t Fold3(const primary_t& r, fold_state_t* fold_state = nullptr);
computed_t Fold2(const primary_t& r, fold_state_t* fold_state = nullptr);
computed_t Fold1(const primary_t& r, fold_state_t* fold_state = nullptr);
computed_t Fold0(const primary_t& r, fold_state_t* fold_state = nullptr);
computed_t FoldBruteForce(const primary_t& r, fold_state_t* fold_state = nullptr);

typedef computed_t (fold_fn_t)(const primary_t&, fold_state_t* fold_state);

fold_fn_t* const FOLD_FUNCTIONS[] = {&Fold0, &Fold1, &Fold2, &Fold3};

const std::map<std::string, ArgParse::option_t> FOLD_OPTIONS = {
    {"alg", ArgParse::option_t("which algorithm for memerna").Arg("0", {"0", "1", "2", "3", "brute"})}
};

fold::fold_fn_t* FoldFunctionFromArgParse(const ArgParse& argparse);

}
}

#endif //MEMERNA_FOLD_H
