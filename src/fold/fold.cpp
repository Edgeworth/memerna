#include <stack>
#include "parsing.h"
#include "fold/fold.h"
#include "fold/traceback.h"
#include "fold/fold_globals.h"
#include "fold/fold_internal.h"

namespace memerna {
namespace fold {

template<typename T>
computed_t FoldInternal(const primary_t& r, fold_state_t* fold_state, T ComputeTables) {
  internal::InitFold(r);
  auto arr = ComputeTables(r);
  traceback_stack_t q;
  auto exterior = TraceExterior(r, arr, q);
  auto computed = TraceStructure(r, arr, exterior, q);
  if (fold_state) {
    fold_state->dp_table = std::move(arr);
    fold_state->ext_table = std::move(exterior);
  }
  return computed;
}

computed_t Fold0(const primary_t& r, fold_state_t* fold_state) {
  return FoldInternal(r, fold_state, internal::ComputeTables0);
}

computed_t Fold1(const primary_t& r, fold_state_t* fold_state) {
  return FoldInternal(r, fold_state, internal::ComputeTables1);
}

computed_t Fold2(const primary_t& r, fold_state_t* fold_state) {
  return FoldInternal(r, fold_state, internal::ComputeTables2);
}

computed_t Fold3(const primary_t& r, fold_state_t* fold_state) {
  return FoldInternal(r, fold_state, internal::ComputeTables3);
}

fold::fold_fn_t* FoldFunctionFromArgParse(const ArgParse& argparse) {
  fold::fold_fn_t* fold_fn = nullptr;
  auto opt = argparse.GetOption("alg");
  if (opt == "0")
    fold_fn = &fold::Fold0;
  else if (opt == "1")
    fold_fn = &fold::Fold1;
  else if (opt == "2")
    fold_fn = &fold::Fold2;
  else if (opt == "3")
    fold_fn = &fold::Fold3;
  else if (opt == "brute")
    fold_fn = &fold::FoldBruteForce;
  else
    verify_expr(false, "unknown fold option");
  return fold_fn;
}

}
}
