// Copyright 2021 Eliot Courtney.
#include "fuzz/fuzzer.h"

#include <cmath>
#include <memory>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

#include "compute/boltz_dp.h"
#include "compute/brute/alg.h"
#include "compute/energy/energy.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "compute/subopt/config.h"
#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "ctx/config.h"
#include "ctx/ctx.h"
#include "model/ctd.h"
#include "model/model.h"
#include "model/secondary.h"
#include "util/array.h"
#include "util/float.h"
#include "util/string.h"

namespace mrna::fuzz {

inline bool equ(BoltzEnergy a, BoltzEnergy b) { return fabs(a - b) < EP; }

Fuzzer::Fuzzer(Primary r, energy::EnergyModelPtr em, FuzzCfg cfg)
    : r_(std::move(r)), em_(em), cfg_(std::move(cfg)) {}

Error Fuzzer::Run() {
  Register("memerna:", MemernaComputeAndCheckState());
  if (cfg_.mfe_rnastructure) Register("rnastructure:", RnastructureComputeAndCheckState());
  if (cfg_.mfe_table) Register("dp tables:", CheckDpTables());
  if (cfg_.subopt) Register("suboptimal:", CheckSuboptimal());

  if (static_cast<int>(r_.size()) <= cfg_.brute_max) Register("brute force:", CheckBruteForce());

  if (cfg_.part) Register("partition:", CheckPartition());

  Register(
      sfmt("Difference on len %zu RNA %s:", r_.size(), r_.ToString().c_str()), std::move(errors_));

  return std::move(errors_);
}

void Fuzzer::Register(const std::string& header, Error&& local) {
  if (local.empty()) return;
  local.push_front(header);
  for (auto& error : local) errors_.push_back("  " + error);
}

bool Fuzzer::HasDuplicates(const std::vector<subopt::SuboptResult>& subopts) {
  struct SuboptCmp {
    bool operator()(const subopt::SuboptResult& a, const subopt::SuboptResult& b) const {
      return std::tie(a.energy, a.tb.s, a.tb.ctd) < std::tie(b.energy, b.tb.s, b.tb.ctd);
    }
  };

  // If energies are different but everything else is the same, it is still a bug.
  // TODO: Sort and avoid copy?
  std::set<subopt::SuboptResult, SuboptCmp> subopt_set;
  for (const auto& subopt : subopts) {
    if (subopt_set.count(subopt)) return true;
    subopt_set.insert(subopt);
  }
  return false;
}

Error Fuzzer::CheckSuboptimalResult(
    const std::vector<subopt::SuboptResult>& subopt, bool has_ctds) {
  Error errors;
  // Check at least one suboptimal structure.
  if (subopt.empty()) errors.push_back("no structures returned");
  // Check MFE.
  if (!subopt.empty() && memerna_subopts_[0].energy != subopt[0].energy)
    errors.push_back(
        sfmt("lowest structure energy %d != mfe %d", subopt[0].energy, memerna_subopts_[0].energy));

  // Only ones with CTDs set can do these tests.
  if (has_ctds) {
    // Check for duplicate structures.
    if (HasDuplicates(subopt)) errors.push_back("has duplicates");

    for (int i = 0; i < static_cast<int>(subopt.size()); ++i) {
      const auto& sub = subopt[i];
      auto suboptimal_efn = em_->TotalEnergy(r_, sub.tb.s, &sub.tb.ctd);
      if (suboptimal_efn.energy != sub.energy) {
        errors.push_back(
            sfmt("structure %d: energy %d != efn %d", i, sub.energy, suboptimal_efn.energy));
        break;
      }

      // Incidentally test ctd parsing.
      auto parsed = ParsePrimaryCtdString(r_.ToString(), sub.tb.ctd.ToString(sub.tb.s));
      if (std::get<Primary>(parsed) != r_ || std::get<Secondary>(parsed) != sub.tb.s ||
          std::get<Ctds>(parsed) != sub.tb.ctd) {
        errors.push_back(sfmt("structure %d: bug in parsing code", i));
        break;
      }
    }
  }
  return errors;
}

Error Fuzzer::CheckSuboptimalResultPair(
    const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b) {
  Error errors;
  if (a.size() != b.size()) {
    errors.push_back(
        sfmt("first has %zu structures != second has %zu structures", a.size(), b.size()));
  } else {
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
      if (a[i].energy != b[i].energy) {
        errors.push_back(sfmt("structure %d: first %d != second %d", i, a[i].energy, b[i].energy));
        break;
      }
    }
  }
  return errors;
}

Error Fuzzer::CheckSuboptimal() {
  Error errors;
  std::vector<std::vector<subopt::SuboptResult>> memerna_subopts_delta, memerna_subopts_max;
  for (auto subopt_alg : ctx::CtxCfg::SUBOPT_ALGS) {
    ctx::CtxCfg cfg{.subopt_alg = subopt_alg};
    ctx::Ctx ctx(em_, cfg);
    memerna_subopts_delta.push_back(
        ctx.SuboptimalIntoVector(Primary(r_), {.delta = cfg_.subopt_delta, .sorted = true}));
    memerna_subopts_max.push_back(
        ctx.SuboptimalIntoVector(Primary(r_), {.strucs = cfg_.subopt_max, .sorted = true}));
  }

  for (int i = 0; i < static_cast<int>(memerna_subopts_delta.size()); ++i) {
    Register(sfmt("memerna delta suboptimal %d:", i),
        CheckSuboptimalResult(memerna_subopts_delta[i], true));
    Register(sfmt("memerna 0 vs memerna %d delta suboptimal:", i),
        CheckSuboptimalResultPair(memerna_subopts_delta[0], memerna_subopts_delta[i]));
  }

  for (int i = 0; i < static_cast<int>(memerna_subopts_max.size()); ++i) {
    Register(
        sfmt("memerna num suboptimal %d:", i), CheckSuboptimalResult(memerna_subopts_max[i], true));
    Register(sfmt("memerna 0 vs memerna %d num suboptimal:", i),
        CheckSuboptimalResultPair(memerna_subopts_max[0], memerna_subopts_max[i]));
  }

#ifdef USE_RNASTRUCTURE
  if (cfg_.subopt_rnastructure) {
    // Suboptimal folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does
    // strange things
    // when the energy for suboptimal structures is 0 or above.
    if (memerna_subopts_[0].energy < -cfg_.subopt_delta) {
      const auto rnastructure_subopt =
          rnastructure_->SuboptimalIntoVector(Primary(r_), cfg_.subopt_delta);
      Register("rnastructure suboptimal:", CheckSuboptimalResult(rnastructure_subopt, false));
      Register("memerna vs rnastructure suboptimal:",
          CheckSuboptimalResultPair(memerna_subopts_delta[0], rnastructure_subopt));
    }
  }
#endif  // USE_RNASTRUCTURE

  return errors;
}

Error Fuzzer::CheckDpTables() {
  Error errors;
  for (int st = static_cast<int>(r_.size()) - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < static_cast<int>(r_.size()); ++en) {
      for (int a = 0; a < DP_SIZE; ++a) {
        const auto memerna0 = memerna_dps[0][st][en][a];
        for (int i = 0; i < static_cast<int>(memerna_dps.size()); ++i) {
          const auto memernai = memerna_dps[i][st][en][a];
          // If meant to be infinity and not.
          if (((memerna0 < CAP_E) != (memernai < CAP_E)) ||
              (memerna0 < CAP_E && memerna0 != memernai)) {
            errors.push_back(sfmt("memerna %d at %d %d %d: %d != %d", i, st, en, a,
                memerna_dps[i][st][en][a], memerna_dps[0][st][en][a]));
            goto loopend;
          }
        }

#ifdef USE_RNASTRUCTURE
        if (cfg_.mfe_rnastructure && (a == DP_P || a == DP_U)) {
          Energy rnastructureval = a == DP_P ? rnastructure_dp_.v.f(st + 1, en + 1)
                                             : rnastructure_dp_.w.f(st + 1, en + 1);
          if (((memerna0 < CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
                  (memerna0 < CAP_E && memerna0 != rnastructureval))) {
            errors.push_back(sfmt("rnastructure at %d %d %d: %d != %d", st, en, a, rnastructureval,
                memerna_dps[0][st][en][a]));
            goto loopend;
          }
        }
#endif  // USE_RNASTRUCTURE
      }
    }
  }
loopend:
  return errors;
}

Error Fuzzer::MemernaComputeAndCheckState() {
  Error errors;
  // Memerna.
  std::vector<Energy> memerna_ctd_efns;
  std::vector<Energy> memerna_optimal_efns;
  for (auto dp_alg : ctx::CtxCfg::DP_ALGS) {
    ctx::Ctx ctx(em_, ctx::CtxCfg{.dp_alg = dp_alg});
    auto res = ctx.Fold(Primary(r_));
    memerna_dps.emplace_back(std::move(res.mfe.dp));
    // First compute with the CTDs that fold returned to check the energy.
    memerna_ctd_efns.push_back(em_->TotalEnergy(r_, res.tb.s, &res.tb.ctd).energy);
    // Also check that the optimal CTD configuration has the same energy.
    // Note that it might not be the same, so we can't do an equality check.
    memerna_optimal_efns.push_back(em_->TotalEnergy(r_, res.tb.s, nullptr).energy);
    memerna_subopts_.push_back(subopt::SuboptResult(std::move(res.tb), res.mfe.energy));
  }

  // Check memerna energies.
  for (int i = 0; i < static_cast<int>(memerna_dps.size()); ++i) {
    if (memerna_subopts_[0].energy != memerna_subopts_[i].energy ||
        memerna_subopts_[0].energy != memerna_ctd_efns[i] ||
        memerna_subopts_[0].energy != memerna_optimal_efns[i])
      errors.push_back(
          sfmt("memerna %d: %d (dp) %d (ctd efn) %d (efn) != mfe %d", i, memerna_subopts_[i].energy,
              memerna_ctd_efns[i], memerna_optimal_efns[i], memerna_subopts_[0].energy));
  }

  return errors;
}

Error Fuzzer::RnastructureComputeAndCheckState() {
  Error errors;
#ifdef USE_RNASTRUCTURE
  // TODO: remove rnastructure_dp_ and use Fold here
  auto fold = rnastructure_->FoldAndDpTable(Primary(r_), &rnastructure_dp_);
  auto efn = rnastructure_->Efn(Primary(r_), Secondary(fold.tb.s));
  // TODO: Test CTDs here as well?
  if (memerna_subopts_[0].energy != fold.mfe.energy || memerna_subopts_[0].energy != efn.energy)
    errors.push_back(sfmt("mfe: rnastructure %d (dp), %d (efn) != mfe %d", fold.mfe.energy,
        efn.energy, memerna_subopts_[0].energy));
#endif  // USE_RNASTRUCTURE
  return errors;
}

Error Fuzzer::CheckBruteForce() {
  Error errors;
  ctx::CtxCfg cfg;
  ctx::Ctx ctx(em_, cfg);

  if (cfg_.subopt) {
    // TODO: move stuff in this function to mfe,subopt,partition.
    auto cfg = subopt::SuboptCfg{.strucs = cfg_.subopt_max, .sorted = true};
    auto brute_subopt = brute::SuboptimalBruteForce(Primary(r_), em_, cfg);
    auto memerna_subopt = ctx.SuboptimalIntoVector(Primary(r_), cfg);

    Register("brute suboptimal:", CheckSuboptimalResult(brute_subopt, true));
    Register("memerna suboptimal:", CheckSuboptimalResult(memerna_subopt, true));
    Register(
        "brute vs memerna suboptimal:", CheckSuboptimalResultPair(brute_subopt, memerna_subopt));
  }

  if (cfg_.part) {
    auto memerna_partition = ctx.Partition(Primary(r_));

    auto brute_partition = brute::PartitionBruteForce(Primary(r_), em_);
    // Types for the partition function are meant to be a bit configurable, so use sstream here.
    if (!equ(brute_partition.part.q, memerna_partition.part.q)) {
      std::stringstream sstream;
      sstream << "q: brute partition " << brute_partition.part.q << " != memerna "
              << memerna_partition.part.q
              << "; difference: " << brute_partition.part.q - memerna_partition.part.q;
      errors.push_back(sstream.str());
    }

    for (int st = 0; st < static_cast<int>(r_.size()); ++st) {
      for (int en = 0; en < static_cast<int>(r_.size()); ++en) {
        if (!equ(brute_partition.part.p[st][en], memerna_partition.part.p[st][en])) {
          std::stringstream sstream;
          sstream << "memerna " << st << " " << en << ": " << memerna_partition.part.p[st][en]
                  << " != brute force " << brute_partition.part.p[st][en] << "; difference: "
                  << brute_partition.part.p[st][en] - memerna_partition.part.p[st][en];
          errors.push_back(sstream.str());
        }
      }
    }
  }
  return errors;
}

Error Fuzzer::CheckPartition() {
  Error errors;
  // TODO: Check DP and ext tables?
  std::vector<part::PartResult> memerna_partitions;
  for (auto part_alg : ctx::CtxCfg::PART_ALGS) {
    ctx::Ctx ctx(em_, ctx::CtxCfg{.part_alg = part_alg});
    memerna_partitions.emplace_back(ctx.Partition(Primary(r_)));
  }

  for (int i = 0; i < static_cast<int>(memerna_partitions.size()); ++i) {
    if (!equ(memerna_partitions[i].part.q, memerna_partitions[0].part.q)) {
      std::stringstream sstream;
      sstream << "q: memerna partition " << i << ": " << memerna_partitions[i].part.q
              << " != " << memerna_partitions[0].part.q
              << "; difference: " << memerna_partitions[i].part.q - memerna_partitions[0].part.q;
      errors.push_back(sstream.str());
    }

    for (int st = 0; st < static_cast<int>(r_.size()); ++st) {
      for (int en = 0; en < static_cast<int>(r_.size()); ++en) {
        if (!equ(memerna_partitions[i].part.p[st][en], memerna_partitions[0].part.p[st][en])) {
          std::stringstream sstream;
          sstream << "memerna " << i << " at " << st << " " << en << ": "
                  << memerna_partitions[i].part.p[st][en]
                  << " != " << memerna_partitions[0].part.p[st][en] << "; difference: "
                  << memerna_partitions[i].part.p[st][en] - memerna_partitions[0].part.p[st][en];
          errors.push_back(sstream.str());
        }
      }
    }

    // TODO: Check probability tables
  }

#ifdef USE_RNASTRUCTURE
  if (cfg_.part_rnastructure) {
    auto rnastructure_part = rnastructure_->Partition(Primary(r_));
    // Types for the partition function are meant to be a bit configurable, so use sstream here.
    if (!equ(rnastructure_part.part.q, memerna_partitions[0].part.q)) {
      std::stringstream sstream;
      sstream << "q: rnastructure partition " << rnastructure_part.part.q << " != memerna "
              << memerna_partitions[0].part.q
              << "; difference: " << rnastructure_part.part.q - memerna_partitions[0].part.q;
      errors.push_back(sstream.str());
    }

    for (int st = 0; st < static_cast<int>(r_.size()); ++st) {
      for (int en = 0; en < static_cast<int>(r_.size()); ++en) {
        if (!equ(rnastructure_part.part.p[st][en], memerna_partitions[0].part.p[st][en])) {
          std::stringstream sstream;
          sstream << "memerna " << st << " " << en << ": " << memerna_partitions[0].part.p[st][en]
                  << " != rnastructure " << rnastructure_part.part.p[st][en] << "; difference: "
                  << rnastructure_part.part.p[st][en] - memerna_partitions[0].part.p[st][en];
          errors.push_back(sstream.str());
        }
      }
    }
  }
#endif  // USE_RNASTRUCTURE

  return errors;
}

}  // namespace mrna::fuzz
