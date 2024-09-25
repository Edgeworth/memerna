// Copyright 2022 Eliot Courtney.
#ifndef FUZZ_FUZZ_INVOCATION_H_
#define FUZZ_FUZZ_INVOCATION_H_
#include <deque>
#include <optional>
#include <string>
#include <vector>

#include "api/ctx/ctx.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "fuzz/fuzz_cfg.h"
#include "model/primary.h"

#ifdef USE_RNASTRUCTURE
#include <memory>
#include <utility>

#include "api/bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

using Error = std::deque<std::string>;

class FuzzInvocation {
 public:
  FuzzInvocation(const Primary& r, std::vector<BackendModelPtr> ms, const FuzzCfg& cfg);

  Error Run();

#ifdef USE_RNASTRUCTURE
  void set_rnastructure(std::shared_ptr<bridge::RNAstructure> rstr) {
    verify(ENERGY_PRECISION == 1, "ENERGY_PRECISION must be 1 for RNAstructure");
    rstr_ = std::move(rstr);
  }
#endif  // USE_RNASTRUCTURE

 private:
  Primary r_;
  std::vector<BackendModelPtr> ms_;
  FuzzCfg cfg_;

  // Store assumed to be correct answers for each problem type:
  std::optional<FoldResult> fold_;
  std::vector<subopt::SuboptResult> subopt_{};
  pfn::PfnResult pfn_{};

  Error errors_;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<bridge::RNAstructure> rstr_;

  Error CheckMfeRNAstructure();
  Error CheckSuboptRNAstructure(subopt::SuboptCfg cfg);
  Error CheckPfnRNAstructure();
#endif  // USE_RNASTRUCTURE

  void Register(const std::string& header, Error&& local);

  void EnsureFoldResult();

  Error CheckMfe();

  Error CheckSubopt();

  static bool SuboptDuplicates(const std::vector<subopt::SuboptResult>& subopts);
  Error CheckSuboptResult(const std::vector<subopt::SuboptResult>& subopt, bool has_ctds = true,
      bool check_duplicates = true);
  static Error CheckSuboptResultPair(subopt::SuboptCfg cfg,
      const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b,
      bool has_ctds = true);

  [[nodiscard]] bool PfnPQEq(flt a, flt b) const;
  [[nodiscard]] bool PfnProbEq(flt a, flt b) const;

  void ComparePfn(
      const PfnTables& got, const PfnTables& want, const std::string& name_got, Error& errors);

  Error CheckPfn();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZ_INVOCATION_H_
