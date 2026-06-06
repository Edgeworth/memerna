// Copyright 2021 Eliot Courtney.
#include "fuzz/fuzz_invocation.h"

#include <fmt/core.h>

#include <algorithm>
#include <iterator>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "api/ctx/algorithm.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy.h"
#include "api/mfe.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/pfn.h"
#include "model/secondary.h"
#include "util/error.h"
#include "util/float.h"
#include "util/log.h"
#include "util/string.h"

namespace mrna::fuzz {

using md::base::DP_P;
using md::base::DP_SIZE;
using md::base::DP_U;
using md::base::DpArrayId;

namespace {

constexpr std::size_t SUBOPT_HARD_CAP = 100000;

[[nodiscard]] std::optional<int> CommonDpSize(const mfe::DpState& state) {
  return std::visit(overloaded{
                        [](const md::base::DpState& base) -> std::optional<int> {
                          return As<int>(base.dp.size());
                        },
                        [](const md::stack::DpState& stack) -> std::optional<int> {
                          return As<int>(stack.base.dp.size());
                        },
                        [](const md::canon::DpState& canon) -> std::optional<int> {
                          return As<int>(canon.dp.size());
                        },
                        [](const std::monostate&) -> std::optional<int> { return std::nullopt; },
                    },
      state);
}

[[nodiscard]] std::optional<Energy> CommonDpIndex(
    const mfe::DpState& state, int st, int en, DpArrayId a) {
  return std::visit(overloaded{
                        [=](const md::base::DpState& base) -> std::optional<Energy> {
                          return base.dp[st][en][a];
                        },
                        [=](const md::stack::DpState& stack) -> std::optional<Energy> {
                          return stack.base.dp[st][en][a];
                        },
                        [=](const md::canon::DpState& canon) -> std::optional<Energy> {
                          switch (a) {
                          case DP_P: return canon.dp[st][en][md::canon::DP_P];
                          case DP_U: return canon.dp[st][en][md::canon::DP_U];
                          case md::base::DP_U2: return canon.dp[st][en][md::canon::DP_U2];
                          default: return std::nullopt;
                          }
                        },
                        [](const std::monostate&) -> std::optional<Energy> { return std::nullopt; },
                    },
      state);
}

[[nodiscard]] std::optional<Energy> CommonExtIndex(const mfe::DpState& state, int st) {
  return std::visit(overloaded{
                        [=](const md::base::DpState& base) -> std::optional<Energy> {
                          return base.ext[st][md::base::EXT];
                        },
                        [=](const md::stack::DpState& stack) -> std::optional<Energy> {
                          return stack.base.ext[st][md::base::EXT];
                        },
                        [=](const md::canon::DpState& canon) -> std::optional<Energy> {
                          return canon.ext[st][md::canon::EXT];
                        },
                        [](const std::monostate&) -> std::optional<Energy> { return std::nullopt; },
                    },
      state);
}

void CompareCommonDpState(const mfe::DpState& got, const mfe::DpState& want,
    const std::string& name_got, const std::string& name_want, Error& errors) {
  const auto got_size = CommonDpSize(got);
  const auto want_size = CommonDpSize(want);
  if (!got_size || !want_size) return;
  if (*got_size == 0 || *want_size == 0) return;  // Brute force doesn't generate tables.
  if (*got_size != *want_size) {
    errors.emplace_back("dp size mismatch:");
    errors.push_back(
        fmt::format("  dp {} vs {}: {} != {}", name_got, name_want, *got_size, *want_size));
    return;
  }

  const int N = *want_size;
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (DpArrayId a : {DP_P, DP_U, md::base::DP_U2}) {
        const auto dp = CommonDpIndex(want, st, en, a);
        const auto dpi = CommonDpIndex(got, st, en, a);
        if (!dp || !dpi) continue;
        if (((*dp < CAP_E) != (*dpi < CAP_E)) || (*dp < CAP_E && *dp != *dpi)) {
          errors.emplace_back("dp mismatch:");
          errors.push_back(fmt::format(
              "  dp {} vs {} at {} {} {}: {} != {}", name_got, name_want, st, en, a, *dpi, *dp));
        }
      }
    }
  }

  for (int st = 0; st < N; ++st) {
    const auto ext = CommonExtIndex(want, st);
    const auto exti = CommonExtIndex(got, st);
    if (!ext || !exti) continue;
    if (((*ext < CAP_E) != (*exti < CAP_E)) || (*ext < CAP_E && *ext != *exti)) {
      errors.emplace_back("ext mismatch:");
      errors.push_back(
          fmt::format("ext {} vs {} at {}: {} != {}", name_got, name_want, st, *exti, *ext));
    }
  }
}

[[nodiscard]] Energy SuboptDeltaForLength(int len) {
  // Number of structures is exponential. Heuristic to estimate delta to target a reasonable number
  // of structures.
  Energy delta = Energy::FromFlt(1.0 / (log(len) + 2.0));
  return delta;
}

[[nodiscard]] bool HasDuplicateTracebacks(const std::vector<subopt::SuboptResult>& subopts) {
  std::set<trace::TraceResult> traceback_set;
  for (const auto& subopt : subopts) {
    if (traceback_set.contains(subopt.tb)) return true;
    traceback_set.insert(trace::TraceResult(subopt.tb));
  }
  return false;
}

#ifdef MRNA_USE_RNASTRUCTURE
[[nodiscard]] bool RNAstructureEnergyCfgSupported(
    const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf) {
  return pf.Empty() && cfg.lonely_pairs == erg::EnergyCfg::LonelyPairs::HEURISTIC &&
      cfg.ctd == erg::EnergyCfg::Ctd::ALL;
}
#endif  // MRNA_USE_RNASTRUCTURE

}  // namespace

FuzzInvocation::FuzzInvocation(const Primary& r, std::vector<BackendModelPtr> ms,
    BackendCfg backend_cfg, erg::PseudofreeCfg pf, const FuzzCfg& fuzz_cfg, bool should_log)
    : r_(r), ms_(std::move(ms)), backend_cfg_(std::move(backend_cfg)), pf_(std::move(pf)),
      cfg_(fuzz_cfg), should_log_(should_log) {
  cfg_.subopt_delta = SuboptDeltaForLength(r.size());
  verify(!ms_.empty(), "must provide at least one energy model to fuzz");
}

Error FuzzInvocation::Run() {
  std::optional<FoldBaseline> fold_baseline;
  if (cfg_.mfe || cfg_.subopt) {
    auto [mfe_errors, baseline] = CheckMfe();
    Register("mfe:", std::move(mfe_errors));
    fold_baseline = std::move(baseline);
  }

  if (cfg_.subopt && fold_baseline.has_value()) Register("subopt:", CheckSubopt(*fold_baseline));

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

const BackendModelPtr* FuzzInvocation::FindSuboptModel(
    const erg::EnergyCfg& energy_cfg, subopt::SuboptCfg subopt_cfg) const {
  for (const auto& m : ms_) {
    auto kind = GetBackendKind(m);
    if (!ResolveMfe(kind, /*alg=*/std::nullopt, backend_cfg_, energy_cfg, pf_, /*log=*/nullptr)
            .has_value())
      continue;
    if (ResolveSubopt(kind, /*alg=*/std::nullopt, backend_cfg_, energy_cfg, pf_, subopt_cfg,
            /*log=*/nullptr)
            .has_value())
      return &m;
  }
  return nullptr;
}

std::tuple<Error, std::optional<FuzzInvocation::FoldBaseline>> FuzzInvocation::CheckMfe() {
  const int N = r_.size();
  Error errors;

  // Run memerna folds.
  std::vector<FoldResult> results;
  std::vector<BackendModelPtr> models;
  std::vector<std::string> tags;
  std::vector<Energy> ctd_efns;  // Efn using returned CTDs.
  std::vector<Energy> opt_efns;  // Efn using optimal CTDs.
  for (const auto& m : ms_) {
    const auto kind = GetBackendKind(m);
    auto maybe_run = [&](MfeAlg mfe_alg) {
      std::string reason;
      if (!MfeAlgIsSupported(kind, mfe_alg, backend_cfg_, cfg_.energy_cfg, pf_, &reason)) {
        if (should_log_) loginfo("mfe NOT fuzzed: {}-{}: {}", kind, mfe_alg, reason);
        return;
      }
      if (should_log_) loginfo("mfe fuzzed: {}-{}", kind, mfe_alg);

      const Ctx ctx(m, backend_cfg_);
      auto res = ctx.Fold(r_, mfe_alg, cfg_.energy_cfg, pf_, {});
      // First compute with the CTDs that fold returned to check the energy.
      ctd_efns.push_back(TotalEnergy(m, r_, res.tb.s, &res.tb.ctd, cfg_.energy_cfg, pf_).energy);

      // Also check that the optimal CTD configuration has the same energy.
      // Note that it might not be the same, so we can't do an peqality check
      // of CTD structure.
      opt_efns.push_back(
          TotalEnergy(m, r_, res.tb.s, /*given_ctd=*/nullptr, cfg_.energy_cfg, pf_).energy);
      results.emplace_back(std::move(res));
      models.push_back(m);
      tags.push_back(fmt::format("{}-{}", kind, mfe_alg));
    };

    for (const auto& entry : MfePriorityForBackend(kind, /*alg=*/std::nullopt, N <= cfg_.brute_max))
      maybe_run(entry.alg);
  }

  if (results.empty()) return {std::move(errors), std::nullopt};

  // Prefer a baseline with a DP table so table comparisons and RNAstructure checks
  // can use the chosen baseline directly.
  int cmp_idx = 0;
  for (int i = 0; i < static_cast<int>(results.size()); ++i) {
    if (mfe::MaybeGetBaseDpState(results[i].mfe.dp) != nullptr) {
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
      errors.push_back(fmt::format("  {}: {} (dp) {} (ctd efn) {} (opt efn) != {} mfe {}", tags[i],
          results[i].mfe.energy, ctd_efns[i], opt_efns[i], tags[cmp_idx], cmp_res.mfe.energy));
    }

    if (cfg_.mfe_table)
      CompareCommonDpState(results[i].mfe.dp, cmp_res.mfe.dp, tags[i], tags[cmp_idx], errors);
  }

  std::optional<FoldBaseline> baseline =
      FoldBaseline{.fold = std::move(cmp_res), .model = std::move(models[cmp_idx])};

#ifdef MRNA_USE_RNASTRUCTURE
  if (cfg_.mfe_rnastructure) {
    Register("RNAstructure:", CheckMfeRNAstructure(baseline.value()));
  }
#endif  // MRNA_USE_RNASTRUCTURE

  return {std::move(errors), std::move(baseline)};
}

Error FuzzInvocation::CheckSubopt(const FoldBaseline& baseline) {
  struct SuboptRun {
    BackendModelPtr model;
    std::vector<subopt::SuboptResult> results;
    std::string tag;
    bool hit_hard_cap = false;
  };

  const int N = r_.size();
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

  for (int cfg_idx = 0; cfg_idx < static_cast<int>(std::size(cfgs)); ++cfg_idx) {
    const auto& cfg = cfgs[cfg_idx];
    auto desc = fmt::format("subopt delta: {} strucs: {} sorted: {}, idx: {}", cfg.delta,
        cfg.strucs, cfg.sorted, cfg_idx);

    auto collect = [&](const BackendModelPtr& m, SuboptAlg subopt_alg) -> std::optional<SuboptRun> {
      const auto kind = GetBackendKind(m);
      std::string reason;
      if (!SuboptAlgIsSupported(
              kind, subopt_alg, backend_cfg_, cfg_.energy_cfg, pf_, cfg, &reason)) {
        if (should_log_) loginfo("subopt NOT fuzzed: {}-{}: {}", kind, subopt_alg, reason);
        return std::nullopt;
      }
      if (should_log_) loginfo("subopt fuzzed: {}-{}", kind, subopt_alg);

      SuboptRun run{
          .model = m, .results = {}, .tag = fmt::format("{}-{}-{}", kind, subopt_alg, cfg_idx)};
      run.results.reserve(SUBOPT_HARD_CAP);
      const Ctx ctx(m, backend_cfg_);
      std::ignore = ctx.Subopt(r_, /*mfe_alg=*/std::nullopt, subopt_alg, cfg_.energy_cfg, pf_, cfg,
          [&](const subopt::SuboptResult& res) {
            run.results.push_back(res);
            run.hit_hard_cap = run.results.size() >= SUBOPT_HARD_CAP;
            return !run.hit_hard_cap;
          });

      // Check sorted if it was meant to be.
      if (cfg.sorted &&
          !std::is_sorted(run.results.begin(), run.results.end(),
              [](const auto& a, const auto& b) { return a.energy < b.energy; }))
        errors.emplace_back("subopt result not sorted by energy");

      std::sort(run.results.begin(), run.results.end());
      return run;
    };

    std::vector<SuboptRun> runs;
    for (const auto& m : ms_) {
      const auto kind = GetBackendKind(m);
      for (const auto& entry :
          SuboptPriorityForBackend(kind, /*alg=*/std::nullopt, /*include_brute=*/false)) {
        // DEBUG algo may OOM even with hard-cap, so only run it (later) if we didn't hit cap.
        if (entry.alg == SuboptAlg::DEBUG) continue;
        if (auto run = collect(m, entry.alg)) runs.push_back(std::move(*run));
      }
    }

    if (runs.empty()) {
      if (should_log_) loginfo("subopt skipped: no phase-1 run supported for {}", desc);
      continue;
    }

    std::vector<std::string> capped_tags;
    std::vector<std::string> completed_tags;
    for (const auto& run : runs)
      (run.hit_hard_cap ? capped_tags : completed_tags).push_back(run.tag);

    if (!capped_tags.empty() && !completed_tags.empty()) {
      errors.emplace_back("subopt hard cap mismatch:");
      errors.push_back(fmt::format("  cfg: {}, capped: {}, completed: {}", desc,
          Join(capped_tags, ","), Join(completed_tags, ",")));
      continue;
    }

    // All hit hard cap, don't run the deferred ones that might OOM.
    if (completed_tags.empty()) continue;

    std::vector<SuboptRun> deferred_runs;
    for (const auto& m : ms_) {
      if (auto run = collect(m, SuboptAlg::DEBUG)) deferred_runs.push_back(std::move(*run));
      if (N <= cfg_.brute_max)
        if (auto run = collect(m, SuboptAlg::BRUTE)) deferred_runs.push_back(std::move(*run));
    }

    bool deferred_hit_cap = false;
    for (auto& run : deferred_runs) {
      if (run.hit_hard_cap) {
        deferred_hit_cap = true;
        errors.emplace_back("subopt deferred hard cap:");
        errors.push_back(fmt::format("  cfg: {}, capped: {}", desc, run.tag));
      }
      runs.push_back(std::move(run));
    }
    if (deferred_hit_cap) continue;

    for (int alg = 0; alg < static_cast<int>(runs.size()); ++alg) {
      Register(fmt::format("{}, cfg: {}", runs[alg].tag, desc),
          CheckSuboptResult(baseline.fold.mfe.energy, runs[alg].results, runs[alg].model));
      Register(fmt::format("{} vs {}, cfg: {}", runs[alg].tag, runs[0].tag, desc),
          CheckSuboptResultPair(cfg, runs[0].results, runs[alg].results));
    }

#ifdef MRNA_USE_RNASTRUCTURE
    if (cfg_idx == 0 && cfg_.subopt_rnastructure) {
      if (RNAstructureEnergyCfgSupported(cfg_.energy_cfg, pf_)) {
        Register("rnastructure:", CheckSuboptRNAstructure(cfg, baseline, runs[0].results));
      } else if (should_log_) {
        loginfo("RNAstructure subopt skipped: unsupported energy or pseudofree configuration");
      }
    }
#endif  // MRNA_USE_RNASTRUCTURE
  }

  return errors;
}

Error FuzzInvocation::CheckSuboptResult(Energy mfe_energy,
    const std::vector<subopt::SuboptResult>& subopt, const BackendModelPtr& m, bool has_ctds,
    bool check_duplicates) {
  Error errors;
  // Check at least one suboptimal structure.
  if (subopt.empty()) errors.emplace_back("no structures returned");
  // Check MFE.
  if (!subopt.empty() && mfe_energy != subopt[0].energy)
    errors.push_back(
        fmt::format("lowest structure energy {} != mfe {}", subopt[0].energy, mfe_energy));

  // Check for duplicate structures.
  if (check_duplicates && HasDuplicateTracebacks(subopt)) errors.emplace_back("has duplicates");

  // Only ones with CTDs set can do these tests.
  // TODO(2): Improve this once we have better CTD option support.
  if (has_ctds) {
    for (int i = 0; i < static_cast<int>(subopt.size()); ++i) {
      const auto& sub = subopt[i];
      auto suboptimal_efn = TotalEnergy(m, r_, sub.tb.s, &sub.tb.ctd, cfg_.energy_cfg, pf_);
      if (suboptimal_efn.energy != sub.energy) {
        errors.push_back(
            fmt::format("structure {}: energy {} != efn {}", i, sub.energy, suboptimal_efn.energy));
        break;
      }

      // Incidentally test ctd parsing.
      auto ctd_string = cfg_.energy_cfg.ToCtdString(sub.tb.s, sub.tb.ctd);
      auto parsed = cfg_.energy_cfg.ParseSeqCtdString(r_.ToSeq(), ctd_string);
      if (std::get<Primary>(parsed) != r_) {
        errors.push_back(fmt::format("structure {}: bug in primary parsing code", i));
        break;
      }
      if (std::get<Secondary>(parsed) != sub.tb.s) {
        errors.push_back(fmt::format("structure {}: bug in secondary parsing code", i));
        break;
      }
      if (std::get<Ctds>(parsed) != sub.tb.ctd) {
        errors.push_back(fmt::format("structure {}: bug in CTD parsing code", i));
        break;
      }
    }
  }
  return errors;
}

Error FuzzInvocation::CheckSuboptResultPair(subopt::SuboptCfg subopt_cfg,
    const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b,
    bool has_ctds) {
  Error errors;
  if (a.size() != b.size()) {
    errors.push_back(
        fmt::format("first has {} structures != second has {} structures", a.size(), b.size()));
  } else {
    for (int64_t i = 0; i < static_cast<int64_t>(a.size()); ++i) {
      // If we were limited by number of structures and we are on the last energy value,
      // different algorithms may not have put the same subset of structures with
      // that energy value into their result, so break.
      if (subopt_cfg.strucs == static_cast<int64_t>(a.size()) && a[i].energy == a.back().energy)
        break;
      if (a[i].energy != b[i].energy)
        errors.push_back(
            fmt::format("structure {}: first {} != second {}", i, a[i].energy, b[i].energy));
      if (a[i].tb.s != b[i].tb.s)
        errors.push_back(fmt::format("structure {}: secondaries differ", i));
      if (has_ctds && a[i].tb.ctd != b[i].tb.ctd)
        errors.push_back(fmt::format("structure {}: ctds differ", i));
    }
  }
  return errors;
}

bool FuzzInvocation::PfnPQEq(flt a, flt b) const {
  return abs_eq(a, b, cfg_.pfn_pq_abs_ep) || rel_eq(a, b, cfg_.pfn_pq_rel_ep);
}

bool FuzzInvocation::PfnProbEq(flt a, flt b) const {
  return abs_eq(a, b, cfg_.pfn_prob_abs_ep) || rel_eq(a, b, cfg_.pfn_prob_rel_ep);
}

void FuzzInvocation::ComparePfn(const PfnTables& got, const PfnTables& want,
    const std::string& name_got, const std::string& name_want, Error& errors) {
  const auto N = As<int>(want.p.size());
  verify(want.p.size() == want.prob.size(), "bug");
  verify(got.p.size() == want.prob.size(), "bug");
  verify(got.p.size() == want.p.size(), "bug");

  if (!PfnPQEq(got.q, want.q)) {
    errors.push_back(fmt::format(
        "{} q: {} != {} q: {}; diff: {}", name_got, got.q, name_want, want.q, got.q - want.q));
  }

  for (int st = 0; st < N; ++st) {
    for (int en = 0; en < N; ++en) {
      if (!PfnPQEq(got.p[st][en], want.p[st][en])) {
        errors.push_back(fmt::format("{} p at [{}, {}]: {} != {} {}; diff: {}", name_got, st, en,
            got.p[st][en], name_want, want.p[st][en], got.p[st][en] - want.p[st][en]));
      }

      if (!PfnProbEq(got.prob[st][en], want.prob[st][en])) {
        errors.push_back(fmt::format("{} prob at [{}, {}]: {} != {} {}; diff: {}", name_got, st, en,
            got.prob[st][en], name_want, want.prob[st][en], got.prob[st][en] - want.prob[st][en]));
      }
    }
  }
}

Error FuzzInvocation::CheckPfn() {
  const int N = r_.size();
  Error errors;
  const erg::EnergyCfgPfn pfn_energy_cfg{cfg_.energy_cfg};
  std::vector<pfn::PfnResult> results;
  std::vector<std::string> tags;
  for (const auto& m : ms_) {
    const auto kind = GetBackendKind(m);
    auto maybe_run = [&](PfnAlg pfn_alg) {
      std::string reason;
      if (!PfnAlgIsSupported(kind, pfn_alg, backend_cfg_, pfn_energy_cfg, pf_, &reason)) {
        if (should_log_) loginfo("pfn NOT fuzzed: {}-{}: {}", kind, pfn_alg, reason);
        return;
      }
      if (should_log_) loginfo("pfn fuzzed: {}-{}", kind, pfn_alg);

      const Ctx ctx(m, backend_cfg_);
      results.emplace_back(ctx.Pfn(r_, pfn_alg, pfn_energy_cfg, pf_));
      tags.push_back(fmt::format("{}-{}", kind, pfn_alg));
    };

    for (const auto& entry : PfnPriorityForBackend(kind, /*alg=*/std::nullopt, N <= cfg_.brute_max))
      maybe_run(entry.alg);
  }

  if (results.empty()) return errors;

  for (int i = 0; i < static_cast<int>(results.size()); ++i)
    ComparePfn(results[i].pfn, results[0].pfn, tags[i], tags[0], errors);

  if (N < cfg_.pfn_subopt) {
    subopt::SuboptCfg subopt_cfg = {.strucs = 100000, .sorted = false};
    if (const auto* m = FindSuboptModel(pfn_energy_cfg, subopt_cfg)) {
      const Ctx ctx(*m, backend_cfg_);
      auto subopts = ctx.SuboptIntoVector(
          r_, /*mfe_alg=*/std::nullopt, /*alg=*/std::nullopt, pfn_energy_cfg, pf_, subopt_cfg);
      flt subopt_q{};
      for (const auto& res : subopts) subopt_q += res.energy.Boltz();

      for (int i = 0; i < static_cast<int>(results.size()); ++i) {
        if (!PfnPQEq(subopt_q, results[i].pfn.q))
          errors.push_back(fmt::format("subopt q: {} != {} pfn q: {}, diff: {}", subopt_q, tags[i],
              results[i].pfn.q, subopt_q - results[i].pfn.q));
      }
    }
  }

  auto pfn_baseline = std::move(results[0]);

#ifdef MRNA_USE_RNASTRUCTURE
  if (cfg_.pfn_rnastructure) Register("RNAstructure:", CheckPfnRNAstructure(pfn_baseline));
#endif  // MRNA_USE_RNASTRUCTURE

  return errors;
}

#ifdef MRNA_USE_RNASTRUCTURE
Error FuzzInvocation::CheckMfeRNAstructure(const FoldBaseline& baseline) {
  const int N = r_.size();
  Error errors;
  dp_state_t rstr_dp;
  auto fold = rstr_->FoldAndDpTable(r_, &rstr_dp);
  auto efn = rstr_->Efn(r_, Secondary(fold.tb.s), erg::EnergyCfg{}, erg::PseudofreeCfg{});

  // Check RNAstructure energies:
  if (baseline.fold.mfe.energy != fold.mfe.energy || baseline.fold.mfe.energy != efn.energy) {
    errors.emplace_back("mfe/efn energy mismatch:");
    errors.push_back(fmt::format(
        "  {} (dp), {} (efn) != mfe {}", fold.mfe.energy, efn.energy, baseline.fold.mfe.energy));
  }

  // Check RNAstructure produced structure:
  // TODO(2): We don't currently pull CTDs from RNAstructure. Also need to
  // rework the efn api to support different CTD options.
  // Also check that the optimal CTD configuration has the same energy.
  // Note that it might not be the same, so we can't do an peqality check
  // of CTD structure.
  auto opt_efn =
      TotalEnergy(baseline.model, r_, fold.tb.s, /*given_ctd=*/nullptr, cfg_.energy_cfg, pf_)
          .energy;
  if (opt_efn != fold.mfe.energy) {
    errors.emplace_back("mfe/efn energy mismatch:");
    errors.push_back(fmt::format("  {} (opt efn) != mfe {}", opt_efn, fold.mfe.energy));
  }

  const auto* want = mfe::MaybeGetBaseDpState(baseline.fold.mfe.dp);
  verify(want != nullptr, "fuzzing with RNAstructure should have base dp state");

  // Check RNAstructure dp table:
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (DpArrayId a = 0; a < DP_SIZE; ++a) {
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

Error FuzzInvocation::CheckSuboptRNAstructure(subopt::SuboptCfg subopt_cfg,
    const FoldBaseline& baseline, const std::vector<subopt::SuboptResult>& subopt) {
  Error errors;
  // Subopt folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does
  // strange things when the energy for suboptimal structures is 0 or above.
  if (subopt.empty()) return errors;
  if (subopt[0].energy < -subopt_cfg.delta) {
    std::vector<subopt::SuboptResult> rstr_subopt;
    rstr_subopt.reserve(SUBOPT_HARD_CAP);

    bool hit_hard_cap = false;
    std::ignore = rstr_->Subopt(r_, /*mfe_alg=*/std::nullopt, /*alg=*/std::nullopt,
        erg::EnergyCfg{}, erg::PseudofreeCfg{}, {.delta = subopt_cfg.delta},
        [&](const subopt::SuboptResult& res) {
          rstr_subopt.push_back(res);
          hit_hard_cap = rstr_subopt.size() >= SUBOPT_HARD_CAP;
          return !hit_hard_cap;
        });

    // If we hit a hard cap, don't compare results, just ignore this run.
    if (hit_hard_cap) return errors;

    std::sort(rstr_subopt.begin(), rstr_subopt.end());
    Register("subopt:",
        CheckSuboptResult(baseline.fold.mfe.energy, rstr_subopt, baseline.model,
            /*has_ctds=*/false, /*check_duplicates=*/false));
    Register("subopt vs memerna:",
        CheckSuboptResultPair(subopt_cfg, subopt, rstr_subopt, /*has_ctds=*/false));
  }

  return errors;
}

Error FuzzInvocation::CheckPfnRNAstructure(const pfn::PfnResult& pfn) {
  Error errors;
  auto rstr_pfn = rstr_->Pfn(r_, /*alg=*/std::nullopt, erg::EnergyCfgPfn{}, erg::PseudofreeCfg{});

  ComparePfn(rstr_pfn.pfn, pfn.pfn, "RNAstructure", "memerna", errors);

  return errors;
}
#endif  // MRNA_USE_RNASTRUCTURE

}  // namespace mrna::fuzz
