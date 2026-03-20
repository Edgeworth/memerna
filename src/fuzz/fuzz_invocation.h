// Copyright 2022 Eliot Courtney.
#ifndef FUZZ_FUZZ_INVOCATION_H_
#define FUZZ_FUZZ_INVOCATION_H_
#include <deque>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "api/ctx/ctx.h"
#include "api/energy/pseudofree_cfg.h"
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
  FuzzInvocation(const Primary& r, std::vector<BackendModelPtr> ms, BackendCfg backend_cfg,
      erg::PseudofreeCfg pf, const FuzzCfg& fuzz_cfg, bool should_log);

  Error Run();

#ifdef USE_RNASTRUCTURE
  void set_rnastructure(std::shared_ptr<bridge::RNAstructure> rstr) {
    verify(ENERGY_PRECISION == 1, "ENERGY_PRECISION must be 1 for RNAstructure");
    rstr_ = std::move(rstr);
  }
#endif  // USE_RNASTRUCTURE

 private:
  struct FoldBaseline {
    FoldResult fold;
    BackendModelPtr model;
  };

  Primary r_;
  std::vector<BackendModelPtr> ms_;
  BackendCfg backend_cfg_;
  erg::PseudofreeCfg pf_;
  FuzzCfg cfg_;
  bool should_log_;

  Error errors_;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<bridge::RNAstructure> rstr_;

  Error CheckMfeRNAstructure(const FoldBaseline& baseline);
  Error CheckSuboptRNAstructure(subopt::SuboptCfg subopt_cfg, const FoldBaseline& baseline,
      const std::vector<subopt::SuboptResult>& subopt);
  Error CheckPfnRNAstructure(const pfn::PfnResult& pfn);
#endif  // USE_RNASTRUCTURE

  void Register(const std::string& header, Error&& local);

  [[nodiscard]] const BackendModelPtr* FindSuboptModel(subopt::SuboptCfg subopt_cfg) const;

  std::tuple<Error, std::optional<FoldBaseline>> CheckMfe();

  Error CheckSubopt(const FoldBaseline& baseline);

  static bool SuboptDuplicates(const std::vector<subopt::SuboptResult>& subopts);
  Error CheckSuboptResult(Energy mfe_energy, const std::vector<subopt::SuboptResult>& subopt,
      const BackendModelPtr& m, bool has_ctds = true, bool check_duplicates = true);
  static Error CheckSuboptResultPair(subopt::SuboptCfg subopt_cfg,
      const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b,
      bool has_ctds = true);

  [[nodiscard]] bool PfnPQEq(flt a, flt b) const;
  [[nodiscard]] bool PfnProbEq(flt a, flt b) const;

  void ComparePfn(const PfnTables& got, const PfnTables& want, const std::string& name_got,
      const std::string& name_want, Error& errors);

  Error CheckPfn();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZ_INVOCATION_H_
