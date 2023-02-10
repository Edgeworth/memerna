// Copyright 2021 Eliot Courtney.
#include "fuzz/fuzz_invocation.h"

#include <fmt/core.h>

#include <algorithm>
#include <compare>
#include <set>
#include <tuple>
#include <utility>
#include <variant>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/energy/energy.h"
#include "api/energy/model.h"
#include "api/mfe.h"
#include "api/part.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/part.h"
#include "model/secondary.h"
#include "models/t04/mfe/dp.h"
#include "util/error.h"
#include "util/float.h"
#include "util/util.h"

namespace mrna::fuzz {

using md::t04::DP_P;
using md::t04::DP_SIZE;
using md::t04::DP_U;
using md::t04::EXT_SIZE;

namespace {

// const flt PROB_EP{0.0001};
// inline bool part_abs_eq(BoltzEnergy a, BoltzEnergy b) { return abs_eq(a, b, PROB_EP); }
inline bool part_rel_eq(BoltzEnergy a, BoltzEnergy b) { return rel_eq(a, b, EP); }

void ComparePart(const Part& got, const Part& want, const std::string& name_got, Error& errors) {
  const int N = static_cast<int>(want.p.size());
  if (!part_rel_eq(got.q, want.q)) {
    errors.push_back(
        fmt::format("{} q: {} != {}; difference: {}", name_got, got.q, want.q, got.q - want.q));
  }

  for (int st = 0; st < N; ++st) {
    for (int en = 0; en < N; ++en) {
      if (!part_rel_eq(got.p[st][en], want.p[st][en])) {
        errors.push_back(fmt::format("{} p at [{}, {}]: {} != {}; difference: {}", name_got, st, en,
            got.p[st][en], want.p[st][en], got.p[st][en] - want.p[st][en]));
      }

      if (!part_rel_eq(got.prob[st][en], want.prob[st][en])) {
        errors.push_back(fmt::format("{} prob at [{}, {}]: {} != {}; difference: {}", name_got, st,
            en, got.prob[st][en], want.prob[st][en], got.prob[st][en] - want.prob[st][en]));
      }
    }
  }
}

void CompareT04DpState(const md::t04::DpState& got, const md::t04::DpState& want,
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
          errors.push_back("dp mismatch:");
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
        errors.push_back("ext mismatch:");
        errors.push_back(fmt::format("ext {} at {} {}: {} != {}", name_got, st, a, exti, ext));
      }
    }
  }
}

}  // namespace

FuzzInvocation::FuzzInvocation(const Primary& r, erg::EnergyModelPtr em, const FuzzCfg& cfg)
    : r_(r), em_(std::move(em)), cfg_(cfg) {}

Error FuzzInvocation::Run() {
  if (cfg_.mfe) Register("mfe:", CheckMfe());
  if (cfg_.subopt) Register("subopt:", CheckSubopt());
  if (cfg_.part) Register("partition:", CheckPartition());

  auto ret = std::move(errors_);
  errors_.clear();
  Register(fmt::format("Difference on len {} RNA {}:", r_.size(), r_.ToSeq()), std::move(ret));

  return std::move(errors_);
}

void FuzzInvocation::Register(const std::string& header, Error&& local) {
  if (local.empty()) return;
  local.push_front(header);
  for (const auto& error : local) errors_.push_back("  " + error);
}

Error FuzzInvocation::CheckMfe() {
  const int N = static_cast<int>(r_.size());
  Error errors;

  // Run memerna folds.
  std::vector<FoldResult> mrna_res;
  std::vector<Energy> mrna_ctd_efns;  // Efn using returned CTDs.
  std::vector<Energy> mrna_opt_efns;  // Efn using optimal CTDs.
  for (auto dp_alg : CtxCfg::DP_ALGS) {
    if (dp_alg == CtxCfg::DpAlg::BRUTE && N > cfg_.brute_max) continue;

    const Ctx ctx(em_, CtxCfg{.dp_alg = dp_alg});
    auto res = ctx.Fold(r_);
    // First compute with the CTDs that fold returned to check the energy.
    mrna_ctd_efns.push_back(erg::TotalEnergy(em_, r_, res.tb.s, &res.tb.ctd).energy);

    // Also check that the optimal CTD configuration has the same energy.
    // Note that it might not be the same, so we can't do an peqality check
    // of CTD structure.
    mrna_opt_efns.push_back(erg::TotalEnergy(em_, r_, res.tb.s, nullptr).energy);
    mrna_res.emplace_back(std::move(res));
  }

  // Check memerna energies compared to themselves and to efn.
  for (int i = 0; i < static_cast<int>(mrna_res.size()); ++i) {
    if (mrna_res[0].mfe.energy != mrna_res[i].mfe.energy ||
        mrna_res[0].mfe.energy != mrna_ctd_efns[i] || mrna_res[0].mfe.energy != mrna_opt_efns[i]) {
      errors.push_back("mfe/efn energy mismatch:");
      errors.push_back(fmt::format("  alg {}: {} (dp) {} (ctd efn) {} (opt efn) != mfe {}", i,
          mrna_res[i].mfe.energy, mrna_ctd_efns[i], mrna_opt_efns[i], mrna_res[0].mfe.energy));
    }

    auto vis = overloaded{
        [&](const md::t04::DpState& got) {
          CompareT04DpState(got, std::get<md::t04::DpState>(mrna_res[0].mfe.dp),
              fmt::format("mrna[{}]", i), errors);
        },
        [&](const md::t22::DpState& got) {
          CompareT04DpState(got.t04, std::get<md::t22::DpState>(mrna_res[0].mfe.dp).t04,
              fmt::format("mrna[{}]", i), errors);
          // TODO(2): test got.penult.
        },
        [&](const std::monostate&) {},
        [&](const auto&) { bug(); },
    };
    std::visit(vis, mrna_res[i].mfe.dp);
  }

  fold_ = std::move(mrna_res[0]);

#ifdef USE_RNASTRUCTURE
  if (cfg_.mfe_rnastructure) Register("rnastructure:", CheckMfeRNAstructure());
#endif  // USE_RNASTRUCTURE

  return errors;
}

Error FuzzInvocation::CheckSubopt() {
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
  std::vector<std::vector<std::vector<subopt::SuboptResult>>> mrna;
  for (auto cfg : cfgs) {
    mrna.emplace_back();
    for (auto subopt_alg : CtxCfg::SUBOPT_ALGS) {
      if (subopt_alg == CtxCfg::SuboptAlg::BRUTE && N > cfg_.brute_max) continue;

      const Ctx ctx(em_, CtxCfg{.subopt_alg = subopt_alg});
      auto res = ctx.SuboptimalIntoVector(r_, cfg);
      // Sort them to make the sorted=false configurations comparable between
      // algoritms.
      std::sort(res.begin(), res.end());
      mrna.back().push_back(std::move(res));
    }
  }

  for (int cfg = 0; cfg < static_cast<int>(mrna.size()); ++cfg) {
    const auto& mrna_cfg = mrna[cfg];
    auto desc = fmt::format("subopt delta: {} strucs: {} sorted: {}", cfgs[cfg].delta,
        cfgs[cfg].strucs, cfgs[cfg].sorted);
    for (int alg = 0; alg < static_cast<int>(mrna_cfg.size()); ++alg) {
      Register(fmt::format("alg {}, cfg: {}", alg, desc), CheckSuboptResult(mrna_cfg[alg], true));
      Register(fmt::format("alg {} vs alg 0, cfg: {}", alg, desc),
          CheckSuboptResultPair(cfgs[cfg], mrna_cfg[0], mrna_cfg[alg]));
    }
  }

  // Put regular configuration (delta-sorted) into common result:
  subopt_ = std::move(mrna.front().front());

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
  Error errors;
  // Check at least one suboptimal structure.
  if (subopt.empty()) errors.push_back("no structures returned");
  // Check MFE.
  if (!subopt.empty() && fold_.mfe.energy != subopt[0].energy)
    errors.push_back(
        fmt::format("lowest structure energy {} != mfe {}", subopt[0].energy, fold_.mfe.energy));

  // Check for duplicate structures.
  if (SuboptDuplicates(subopt)) errors.push_back("has duplicates");

  // Only ones with CTDs set can do these tests.
  // TODO(2): Improve this once we have better CTD option support.
  if (has_ctds) {
    for (int i = 0; i < static_cast<int>(subopt.size()); ++i) {
      const auto& sub = subopt[i];
      auto suboptimal_efn = erg::TotalEnergy(em_, r_, sub.tb.s, &sub.tb.ctd);
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

Error FuzzInvocation::CheckPartition() {
  const int N = static_cast<int>(r_.size());
  Error errors;
  std::vector<part::PartResult> mrna_parts;
  for (auto part_alg : CtxCfg::PART_ALGS) {
    if (part_alg == CtxCfg::PartAlg::BRUTE && N > cfg_.brute_max) continue;

    const Ctx ctx(em_, CtxCfg{.part_alg = part_alg});
    mrna_parts.emplace_back(ctx.Partition(r_));
  }

  for (int i = 0; i < static_cast<int>(mrna_parts.size()); ++i)
    ComparePart(mrna_parts[i].part, mrna_parts[0].part, fmt::format("memerna[{}]", i), errors);

  part_ = std::move(mrna_parts[0]);

#ifdef USE_RNASTRUCTURE
  if (cfg_.part_rnastructure) Register("rnastructure:", CheckPartitionRNAstructure());
#endif  // USE_RNASTRUCTURE

  return errors;
}

#ifdef USE_RNASTRUCTURE
Error FuzzInvocation::CheckMfeRNAstructure() {
  const int N = static_cast<int>(r_.size());
  Error errors;
  dp_state_t rstr_dp;
  auto fold = rstr_->FoldAndDpTable(r_, &rstr_dp);
  auto efn = rstr_->Efn(r_, Secondary(fold.tb.s));

  // Check RNAstructure energies:
  if (fold_.mfe.energy != fold.mfe.energy || fold_.mfe.energy != efn.energy) {
    errors.push_back("mfe/efn energy mismatch:");
    errors.push_back(fmt::format(
        "  {} (dp), {} (efn) != mfe {}", fold.mfe.energy, efn.energy, fold_.mfe.energy));
  }

  // Check RNAstructure produced structure:
  // TODO(2): We don't currently pull CTDs from RNAstructure. Also need to
  // rework the efn api to support different CTD options.
  // Also check that the optimal CTD configuration has the same energy.
  // Note that it might not be the same, so we can't do an peqality check
  // of CTD structure.
  auto opt_efn = erg::TotalEnergy(em_, r_, fold.tb.s, nullptr).energy;
  if (opt_efn != fold.mfe.energy) {
    errors.push_back("mfe/efn energy mismatch:");
    errors.push_back(fmt::format("  {} (opt efn) != mfe {}", opt_efn, fold.mfe.energy));
  }

  auto vis = overloaded{
      [&](const md::t04::DpState& want) {
        // Check RNAstructure dp table:
        for (int st = N - 1; st >= 0; --st) {
          for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
            for (int a = 0; a < DP_SIZE; ++a) {
              auto dp = want.dp[st][en][a];
              if (a == DP_P || a == DP_U) {
                auto rstr_eval =
                    a == DP_P ? rstr_dp.v.f(st + 1, en + 1) : rstr_dp.w.f(st + 1, en + 1);
                if (((dp < CAP_E) != (rstr_eval < INFINITE_ENERGY - 1000) ||
                        (dp < CAP_E && dp != bridge::RNAstructure::ToEnergy(rstr_eval)))) {
                  errors.push_back("dp mismatch:");
                  errors.push_back(
                      fmt::format("  dp at {} {} {}: {} != {}", st, en, a, rstr_eval, dp));
                }
              }
            }
          }
        }
      },
      [&](const auto&) { bug(); },
  };
  std::visit(vis, fold_.mfe.dp);

  // TODO(2): Check RNAstructure ext table.
  return errors;
}

Error FuzzInvocation::CheckSuboptRNAstructure(subopt::SuboptCfg cfg) {
  Error errors;
  // Suboptimal folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does
  // strange things when the energy for suboptimal structures is 0 or above.
  if (subopt_[0].energy < -cfg_.subopt_delta) {
    const auto rstr_subopt = rstr_->SuboptimalIntoVector(r_, cfg_.subopt_delta);
    Register("subopt:", CheckSuboptResult(rstr_subopt, false));
    Register("subopt vs memerna:", CheckSuboptResultPair(cfg, subopt_, rstr_subopt));
  }

  return errors;
}

Error FuzzInvocation::CheckPartitionRNAstructure() {
  Error errors;
  auto rstr_part = rstr_->Partition(r_);

  ComparePart(rstr_part.part, part_.part, "rnastructure", errors);

  return errors;
}
#endif  // USE_RNASTRUCTURE

}  // namespace mrna::fuzz
