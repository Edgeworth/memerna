// Copyright 2021 Eliot Courtney.
#include "fuzz/fuzz_invocation.h"

#include <fmt/core.h>

#include <algorithm>
#include <set>
#include <utility>
#include <variant>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/energy/energy.h"
#include "api/mfe.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "backends/common/base/dp.h"
#include "backends/stack/mfe/mfe.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/pfn.h"
#include "model/secondary.h"
#include "util/error.h"
#include "util/float.h"
#include "util/util.h"

namespace mrna::fuzz {

using md::base::DP_P;
using md::base::DP_SIZE;
using md::base::DP_U;
using md::base::EXT_SIZE;

namespace {

void CompareBaseDpState(const md::base::DpState& got, const md::base::DpState& want,
    const std::string& name_got, Error& errors) {
  if (got.dp.empty()) return;  // Brute force doesn't generate tables.

  const int N = static_cast<int>(want.dp.size());
  // Check dp tables:
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (int a = 0; a < DP_SIZE; ++a) {
        auto dp = want.dp[st][en][a];

        auto dpi = got.dp[st][en][a];
        // If meant to be infinity and not.
        if (((dp < CAP_E) != (dpi < CAP_E)) || (dp < CAP_E && dp != dpi)) {
          errors.emplace_back("dp mismatch:");
          errors.push_back(
              fmt::format("  dp {} at {} {} {}: {} != {}", name_got, st, en, a, dpi, dp));
        }
      }
    }
  }

  // Check ext tables:
  for (int st = 0; st < N; ++st) {
    for (int a = 0; a < EXT_SIZE; ++a) {
      auto ext = want.ext[st][a];
      auto exti = got.ext[st][a];
      // If meant to be infinity and not.
      if (((ext < CAP_E) != (exti < CAP_E)) || (ext < CAP_E && ext != exti)) {
        errors.emplace_back("ext mismatch:");
        errors.push_back(fmt::format("ext {} at {} {}: {} != {}", name_got, st, a, exti, ext));
      }
    }
  }
}

md::base::DpState* MaybeGetBaseDpState(mfe::DpState& dp) {
  auto base_dp = overloaded{
      [&](md::base::DpState& got) -> md::base::DpState* { return &got; },
      [&](md::stack::DpState& got) -> md::base::DpState* { return &got.base; },
      [&](const std::monostate&) -> md::base::DpState* { return nullptr; },
      [&](const auto&) -> md::base::DpState* {
        fatal("bug");
        return nullptr;
      },
  };
  return std::visit(base_dp, dp);
}

}  // namespace

FuzzInvocation::FuzzInvocation(
    const Primary& r, std::vector<BackendModelPtr> ms, const FuzzCfg& cfg)
    : r_(r), ms_(std::move(ms)), cfg_(cfg) {
  verify(!ms_.empty(), "must provide at least one energy model to fuzz");
}

Error FuzzInvocation::Run() {
  if (cfg_.mfe) Register("mfe:", CheckMfe());
  if (cfg_.subopt) Register("subopt:", CheckSubopt());
  if (cfg_.pfn) Register("pfn:", CheckPfn());

  auto ret = std::move(errors_);
  errors_.clear();
  Register(fmt::format("Diff on len {} RNA {}:", r_.size(), r_.ToSeq()), std::move(ret));

  return std::move(errors_);
}

void FuzzInvocation::Register(const std::string& header, Error&& local) {
  auto errors = std::move(local);
  if (errors.empty()) return;
  errors.push_front(header);
  for (const auto& error : errors) errors_.push_back("  " + error);
}

void FuzzInvocation::EnsureFoldResult() {
  if (!fold_) fold_ = Ctx(ms_[0], CtxCfg{}).Fold(r_, {});
}

Error FuzzInvocation::CheckMfe() {
  const int N = static_cast<int>(r_.size());
  Error errors;

  // Run memerna folds.
  std::vector<FoldResult> results;
  std::vector<Energy> ctd_efns;  // Efn using returned CTDs.
  std::vector<Energy> opt_efns;  // Efn using optimal CTDs.
  for (const auto& m : ms_) {
    for (auto mfe_alg : CtxCfg::MfeAlgsForBackend(m)) {
      if (mfe_alg == CtxCfg::MfeAlg::AUTO) continue;
      if (mfe_alg == CtxCfg::MfeAlg::BRUTE && N > cfg_.brute_max) continue;

      const Ctx ctx(m, CtxCfg{.mfe_alg = mfe_alg});
      auto res = ctx.Fold(r_, {});
      // First compute with the CTDs that fold returned to check the energy.
      ctd_efns.push_back(TotalEnergy(m, r_, res.tb.s, &res.tb.ctd).energy);

      // Also check that the optimal CTD configuration has the same energy.
      // Note that it might not be the same, so we can't do an peqality check
      // of CTD structure.
      opt_efns.push_back(TotalEnergy(m, r_, res.tb.s, nullptr).energy);
      results.emplace_back(std::move(res));
    }
  }

  // Find first dp table that exists.
  md::base::DpState* base_dp = nullptr;
  int cmp_idx = 0;
  for (int i = 0; i < static_cast<int>(results.size()); ++i) {
    base_dp = MaybeGetBaseDpState(results[i].mfe.dp);
    if (base_dp) {
      cmp_idx = i;
      break;
    }
  }

  // Check memerna energies compared to themselves and to efn.
  auto& cmp_res = results[cmp_idx];
  for (int i = 0; i < static_cast<int>(results.size()); ++i) {
    if (cmp_res.mfe.energy != results[i].mfe.energy || cmp_res.mfe.energy != ctd_efns[i] ||
        cmp_res.mfe.energy != opt_efns[i]) {
      errors.emplace_back("mfe/efn energy mismatch:");
      errors.push_back(fmt::format("  alg {}: {} (dp) {} (ctd efn) {} (opt efn) != alg {} mfe {}",
          i, results[i].mfe.energy, ctd_efns[i], opt_efns[i], cmp_idx, cmp_res.mfe.energy));
    }

    if (cfg_.mfe_table && base_dp) {
      auto* got = MaybeGetBaseDpState(results[i].mfe.dp);
      if (got) CompareBaseDpState(*got, *base_dp, fmt::format("mrna[{}]", i), errors);
    }
  }

  fold_ = std::move(cmp_res);

#ifdef USE_RNASTRUCTURE
  if (cfg_.mfe_rnastructure) Register("RNAstructure:", CheckMfeRNAstructure());
#endif  // USE_RNASTRUCTURE

  return errors;
}

Error FuzzInvocation::CheckSubopt() {
  EnsureFoldResult();

  const int N = static_cast<int>(r_.size());
  Error errors;

  subopt::SuboptCfg cfgs[] = {
      // Delta-sorted is the default config to compare other packages with.
      {.delta = cfg_.subopt_delta, .sorted = true},
      {.strucs = cfg_.subopt_strucs, .sorted = true},
      {.delta = cfg_.subopt_delta, .sorted = false},
      {.strucs = cfg_.subopt_strucs, .sorted = false},
      {.delta = cfg_.subopt_delta, .strucs = cfg_.subopt_strucs, .sorted = true},
      {.delta = cfg_.subopt_delta, .strucs = cfg_.subopt_strucs, .sorted = false},
  };
  std::vector<std::pair<subopt::SuboptCfg, std::vector<std::vector<subopt::SuboptResult>>>> results;
  for (auto cfg : cfgs) {
    results.push_back({cfg, {}});
    for (const auto& m : ms_) {
      for (auto subopt_alg : CtxCfg::SuboptAlgsForBackend(m)) {
        if (subopt_alg == CtxCfg::SuboptAlg::AUTO) continue;
        if (subopt_alg == CtxCfg::SuboptAlg::BRUTE && N > cfg_.brute_max) continue;

        const Ctx ctx(m, CtxCfg{.subopt_alg = subopt_alg});
        auto res = ctx.SuboptIntoVector(r_, cfg);
        // Sort them to make the sorted=false configurations comparable between
        // algoritms.
        std::sort(res.begin(), res.end());
        results.back().second.push_back(std::move(res));
      }
    }
  }

  for (int i = 0; i < static_cast<int>(results.size()); ++i) {
    const auto& [cfg, res] = results[i];
    auto desc = fmt::format(
        "subopt delta: {} strucs: {} sorted: {}, idx: {}", cfg.delta, cfg.strucs, cfg.sorted, i);
    for (int alg = 0; alg < static_cast<int>(res.size()); ++alg) {
      Register(fmt::format("alg {}, cfg: {}", alg, desc),
          CheckSuboptResult(res[alg], /*has_ctds=*/true));
      Register(fmt::format("alg {} vs alg 0, cfg: {}", alg, desc),
          CheckSuboptResultPair(cfg, res[0], res[alg]));
    }
  }

  // Put regular configuration (delta-sorted) into common result:
  subopt_ = std::move(results.front().second.front());

#ifdef USE_RNASTRUCTURE
  if (cfg_.subopt_rnastructure) Register("rnastructure:", CheckSuboptRNAstructure(cfgs[0]));
#endif  // USE_RNASTRUCTURE

  return errors;
}

bool FuzzInvocation::SuboptDuplicates(const std::vector<subopt::SuboptResult>& subopts) {
  // If energies are different but everything else is the same, it is still a bug.
  std::set<subopt::SuboptResult> subopt_set;
  for (const auto& subopt : subopts) {
    if (subopt_set.contains(subopt)) return true;
    subopt_set.insert(subopt);
  }
  return false;
}

Error FuzzInvocation::CheckSuboptResult(
    const std::vector<subopt::SuboptResult>& subopt, bool has_ctds) {
  verify(fold_.has_value(), "bug");
  Error errors;
  // Check at least one suboptimal structure.
  if (subopt.empty()) errors.emplace_back("no structures returned");
  // Check MFE.
  if (!subopt.empty() && fold_->mfe.energy != subopt[0].energy)
    errors.push_back(
        fmt::format("lowest structure energy {} != mfe {}", subopt[0].energy, fold_->mfe.energy));

  // Check for duplicate structures.
  if (SuboptDuplicates(subopt)) errors.emplace_back("has duplicates");

  // Only ones with CTDs set can do these tests.
  // TODO(2): Improve this once we have better CTD option support.
  if (has_ctds) {
    for (int i = 0; i < static_cast<int>(subopt.size()); ++i) {
      const auto& sub = subopt[i];
      auto suboptimal_efn = TotalEnergy(ms_[0], r_, sub.tb.s, &sub.tb.ctd);
      if (suboptimal_efn.energy != sub.energy) {
        errors.push_back(
            fmt::format("structure {}: energy {} != efn {}", i, sub.energy, suboptimal_efn.energy));
        break;
      }

      // Incidentally test ctd parsing.
      auto parsed = ParseSeqCtdString(r_.ToSeq(), sub.tb.ctd.ToString(sub.tb.s));
      if (std::get<Primary>(parsed) != r_ || std::get<Secondary>(parsed) != sub.tb.s ||
          std::get<Ctds>(parsed) != sub.tb.ctd) {
        errors.push_back(fmt::format("structure {}: bug in parsing code", i));
        break;
      }
    }
  }
  return errors;
}

Error FuzzInvocation::CheckSuboptResultPair(subopt::SuboptCfg cfg,
    const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b) {
  Error errors;
  if (a.size() != b.size()) {
    errors.push_back(
        fmt::format("first has {} structures != second has {} structures", a.size(), b.size()));
  } else {
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
      // If we were limited by number of structures and we are on the last energy value,
      // different algorithms may not have put the same subset of structures with
      // that energy value into their result, so break.
      if (cfg.strucs == static_cast<int>(a.size()) && a[i].energy == a.back().energy) break;
      if (a[i].energy != b[i].energy)
        errors.push_back(
            fmt::format("structure {}: first {} != second {}", i, a[i].energy, b[i].energy));
      if (a[i].tb.s != b[i].tb.s)
        errors.push_back(fmt::format("structure {}: secondaries differ", i));
      if (a[i].tb.ctd != b[i].tb.ctd)
        errors.push_back(fmt::format("structure {}: ctds differ", i, a[i].energy, b[i].energy));
    }
  }
  return errors;
}

bool FuzzInvocation::PfnPQEq(flt a, flt b) const { return absrel_eq(a, b, cfg_.pfn_pq_ep); }

bool FuzzInvocation::PfnProbEq(flt a, flt b) const { return absrel_eq(a, b, cfg_.pfn_prob_ep); }

void FuzzInvocation::ComparePfn(
    const PfnTables& got, const PfnTables& want, const std::string& name_got, Error& errors) {
  const int N = static_cast<int>(want.p.size());
  verify(want.prob.size() == want.prob.size(), "bug");
  verify(got.p.size() == want.prob.size(), "bug");

  if (!PfnPQEq(got.q, want.q)) {
    errors.push_back(
        fmt::format("{} q: {} != {}; diff: {}", name_got, got.q, want.q, got.q - want.q));
  }

  for (int st = 0; st < N; ++st) {
    for (int en = 0; en < N; ++en) {
      if (!PfnPQEq(got.p[st][en], want.p[st][en])) {
        errors.push_back(fmt::format("{} p at [{}, {}]: {} != {}; diff: {}", name_got, st, en,
            got.p[st][en], want.p[st][en], got.p[st][en] - want.p[st][en]));
      }

      if (!PfnProbEq(got.prob[st][en], want.prob[st][en])) {
        errors.push_back(fmt::format("{} prob at [{}, {}]: {} != {}; diff: {}", name_got, st, en,
            got.prob[st][en], want.prob[st][en], got.prob[st][en] - want.prob[st][en]));
      }
    }
  }
}

Error FuzzInvocation::CheckPfn() {
  const int N = static_cast<int>(r_.size());
  Error errors;
  std::vector<pfn::PfnResult> results;
  for (const auto& m : ms_) {
    for (auto pfn_alg : CtxCfg::PfnAlgsForBackend(m)) {
      if (pfn_alg == CtxCfg::PfnAlg::AUTO) continue;
      if (pfn_alg == CtxCfg::PfnAlg::BRUTE && N > cfg_.brute_max) continue;

      const Ctx ctx(m, CtxCfg{.pfn_alg = pfn_alg});
      results.emplace_back(ctx.Pfn(r_));
    }
  }

  for (int i = 0; i < static_cast<int>(results.size()); ++i)
    ComparePfn(results[i].pfn, results[0].pfn, fmt::format("memerna[{}]", i), errors);

  pfn_ = std::move(results[0]);

#ifdef USE_RNASTRUCTURE
  if (cfg_.pfn_rnastructure) Register("RNAstructure:", CheckPfnRNAstructure());
#endif  // USE_RNASTRUCTURE

  return errors;
}

#ifdef USE_RNASTRUCTURE
Error FuzzInvocation::CheckMfeRNAstructure() {
  verify(fold_.has_value(), "bug");

  const int N = static_cast<int>(r_.size());
  Error errors;
  dp_state_t rstr_dp;
  auto fold = rstr_->FoldAndDpTable(r_, &rstr_dp);
  auto efn = rstr_->Efn(r_, Secondary(fold.tb.s));

  // Check RNAstructure energies:
  if (fold_->mfe.energy != fold.mfe.energy || fold_->mfe.energy != efn.energy) {
    errors.emplace_back("mfe/efn energy mismatch:");
    errors.push_back(fmt::format(
        "  {} (dp), {} (efn) != mfe {}", fold.mfe.energy, efn.energy, fold_->mfe.energy));
  }

  // Check RNAstructure produced structure:
  // TODO(2): We don't currently pull CTDs from RNAstructure. Also need to
  // rework the efn api to support different CTD options.
  // Also check that the optimal CTD configuration has the same energy.
  // Note that it might not be the same, so we can't do an peqality check
  // of CTD structure.
  auto opt_efn = TotalEnergy(ms_[0], r_, fold.tb.s, nullptr).energy;
  if (opt_efn != fold.mfe.energy) {
    errors.emplace_back("mfe/efn energy mismatch:");
    errors.push_back(fmt::format("  {} (opt efn) != mfe {}", opt_efn, fold.mfe.energy));
  }

  auto* want = MaybeGetBaseDpState(fold_->mfe.dp);
  verify(want != nullptr, "fuzzing with RNAstructure should have base dp state");

  // Check RNAstructure dp table:
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (int a = 0; a < DP_SIZE; ++a) {
        auto dp = want->dp[st][en][a];
        if (a == DP_P || a == DP_U) {
          auto rstr_eval = a == DP_P ? rstr_dp.v.f(st + 1, en + 1) : rstr_dp.w.f(st + 1, en + 1);
          if (((dp < CAP_E) != (rstr_eval < INFINITE_ENERGY - 1000) ||
                  (dp < CAP_E && dp != bridge::RNAstructure::ToEnergy(rstr_eval)))) {
            errors.emplace_back("dp mismatch:");
            errors.push_back(fmt::format("  dp at {} {} {}: {} != {}", st, en, a, rstr_eval, dp));
          }
        }
      }
    }
  }

  // TODO(2): Check RNAstructure ext table.
  return errors;
}

Error FuzzInvocation::CheckSuboptRNAstructure(subopt::SuboptCfg cfg) {
  Error errors;
  // Subopt folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does
  // strange things when the energy for suboptimal structures is 0 or above.
  if (subopt_[0].energy < -cfg_.subopt_delta) {
    const auto rstr_subopt = rstr_->SuboptIntoVector(r_, cfg_.subopt_delta);
    Register("subopt:", CheckSuboptResult(rstr_subopt, false));
    Register("subopt vs memerna:", CheckSuboptResultPair(cfg, subopt_, rstr_subopt));
  }

  return errors;
}

Error FuzzInvocation::CheckPfnRNAstructure() {
  Error errors;
  auto rstr_pfn = rstr_->Pfn(r_);

  ComparePfn(rstr_pfn.pfn, pfn_.pfn, "RNAstructure", errors);

  return errors;
}
#endif  // USE_RNASTRUCTURE

}  // namespace mrna::fuzz
