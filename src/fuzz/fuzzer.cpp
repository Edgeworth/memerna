// Copyright 2021 E.
#include "fuzz/fuzzer.h"

#include <cinttypes>
#include <set>
#include <sstream>

#include "compute/energy/load_model.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/brute.h"
#include "compute/subopt/brute.h"
#include "model/context.h"
#include "util/string.h"

namespace mrna::fuzz {

using mfe::internal::gdp;
using partition::PartitionBruteForce;
using subopt::SuboptimalBruteForce;

inline bool equ(BoltzEnergy a, BoltzEnergy b) { return fabs(a - b) < EP; }

Fuzzer::Fuzzer(Primary r, FuzzCfg cfg, energy::EnergyModel em)
    : r_(std::move(r)), cfg_(std::move(cfg)),
      em_(cfg_.random_model ? energy::LoadRandomEnergyModel(cfg_.seed) : em) {}

Error Fuzzer::Run() {
  Error errors;
  AppendErrors(errors, MaybePrepend(MemernaComputeAndCheckState(), "memerna:"));
  if (cfg_.mfe_rnastructure)
    AppendErrors(errors, MaybePrepend(RnastructureComputeAndCheckState(), "rnastructure:"));
  if (cfg_.table_check) AppendErrors(errors, MaybePrepend(CheckDpTables(), "dp tables:"));
  if (cfg_.subopt) AppendErrors(errors, MaybePrepend(CheckSuboptimal(), "suboptimal:"));

  if (static_cast<int>(r_.size()) <= cfg_.brute_cutoff)
    AppendErrors(errors, MaybePrepend(CheckBruteForce(), "brute force:"));

  if (cfg_.partition) AppendErrors(errors, MaybePrepend(CheckPartition(), "partition:"));

  if (!errors.empty()) {
    if (cfg_.random_model)
      errors.push_front(sfmt("Used random energy model with seed: %" PRIuFAST32 "\n", cfg_.seed));
    else
      errors.push_front(sfmt("Used specified energy model"));
    errors = MaybePrepend(
        errors, sfmt("Difference on len %zu RNA %s:", r_.size(), PrimaryToString(r_).c_str()));
  }

  return errors;
}

Error Fuzzer::MaybePrepend(const Error& main, const std::string& header) {
  if (main.empty()) return main;
  Error nmain;
  nmain.push_front(header);
  for (auto& error : main) nmain.push_back("  " + error);  // mfw this inefficiency
  return nmain;
}

void Fuzzer::AppendErrors(Error& main, Error&& extra) {
  for (auto& s : extra) main.push_back(std::move(s));
}

bool Fuzzer::HasDuplicates(const std::vector<Computed>& computeds) {
  // If energies are different but everything else is the same, it is still a bug.
  std::set<std::pair<Secondary, std::vector<Ctd>>> suboptimal_set;
  for (const auto& computed : computeds) {
    auto val = std::make_pair(computed.s, computed.base_ctds);
    if (suboptimal_set.count(val)) return true;
    suboptimal_set.insert(val);
  }
  return false;
}

Error Fuzzer::CheckSuboptimalResult(const std::vector<Computed>& subopt, bool has_ctds) {
  Error errors;
  // Check at least one suboptimal structure.
  if (subopt.empty()) errors.push_back("no structures returned");
  // Check MFE.
  if (!subopt.empty() && memerna_computeds_[0].energy != subopt[0].energy)
    errors.push_back(sfmt(
        "lowest structure energy %d != mfe %d", subopt[0].energy, memerna_computeds_[0].energy));

  // Only ones with CTDs set can do these tests.
  if (has_ctds) {
    // Check for duplicate structures.
    if (HasDuplicates(subopt)) errors.push_back("has duplicates");

    for (int i = 0; i < static_cast<int>(subopt.size()); ++i) {
      const auto& structure = subopt[i];
      auto suboptimal_efn = energy::ComputeEnergyWithCtds(structure, em_);
      if (suboptimal_efn.energy != structure.energy) {
        errors.push_back(
            sfmt("structure %d: energy %d != efn %d", i, structure.energy, suboptimal_efn.energy));
        break;
      }

      // Incidentally test ctd parsing.
      auto parsed_computed =
          ParseCtdComputed(PrimaryToString(structure.s.r), ComputedToCtdString(structure));
      parsed_computed.energy = structure.energy;
      if (parsed_computed != structure) {
        errors.push_back(sfmt("structure %d: bug in parsing code", i));
        break;
      }
    }
  }
  return errors;
}

Error Fuzzer::CheckSuboptimalResultPair(
    const std::vector<Computed>& a, const std::vector<Computed>& b) {
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
  std::vector<std::vector<Computed>> memerna_subopts_delta, memerna_subopts_num;
  for (auto subopt_alg : ModelCfg::SUBOPTIMAL_ALGS) {
    ModelCfg cfg(ModelCfg::TableAlg::TWO, subopt_alg);
    Context ctx(r_, em_, cfg);
    memerna_subopts_delta.push_back(ctx.SuboptimalIntoVector(true, cfg_.subopt_delta, -1));
    memerna_subopts_num.push_back(ctx.SuboptimalIntoVector(true, -1, cfg_.subopt_max));
  }

  for (int i = 0; i < static_cast<int>(memerna_subopts_delta.size()); ++i) {
    AppendErrors(errors,
        MaybePrepend(CheckSuboptimalResult(memerna_subopts_delta[i], true),
            sfmt("memerna delta suboptimal %d:", i)));
    AppendErrors(errors,
        MaybePrepend(CheckSuboptimalResultPair(memerna_subopts_delta[0], memerna_subopts_delta[i]),
            sfmt("memerna 0 vs memerna %d delta suboptimal:", i)));
  }

  for (int i = 0; i < static_cast<int>(memerna_subopts_num.size()); ++i) {
    AppendErrors(errors,
        MaybePrepend(CheckSuboptimalResult(memerna_subopts_num[i], true),
            sfmt("memerna num suboptimal %d:", i)));
    AppendErrors(errors,
        MaybePrepend(CheckSuboptimalResultPair(memerna_subopts_num[0], memerna_subopts_num[i]),
            sfmt("memerna 0 vs memerna %d num suboptimal:", i)));
  }

#ifdef USE_RNASTRUCTURE
  if (cfg_.subopt_rnastructure) {
    // Suboptimal folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does
    // strange things
    // when the energy for suboptimal structures is 0 or above.
    if (memerna_computeds_[0].energy < -cfg_.subopt_delta) {
      const auto rnastructure_subopt = rnastructure_->SuboptimalIntoVector(r_, cfg_.subopt_delta);
      AppendErrors(errors,
          MaybePrepend(
              CheckSuboptimalResult(rnastructure_subopt, false), "rnastructure suboptimal:"));
      AppendErrors(errors,
          MaybePrepend(CheckSuboptimalResultPair(memerna_subopts_delta[0], rnastructure_subopt),
              "memerna vs rnastructure suboptimal:"));
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
  for (auto table_alg : ModelCfg::TABLE_ALGS) {
    Context ctx(r_, em_, ModelCfg(table_alg));
    auto computed = ctx.Fold();
    memerna_dps.emplace_back(std::move(gdp));
    // First compute with the CTDs that fold returned to check the energy.
    memerna_ctd_efns.push_back(energy::ComputeEnergyWithCtds(computed, em_).energy);
    // Also check that the optimal CTD configuration has the same energy.
    // Note that it might not be the same, so we can't do an equality check.
    memerna_optimal_efns.push_back(energy::ComputeEnergy(computed.s, em_).energy);
    memerna_computeds_.push_back(std::move(computed));
  }

  // Check memerna energies.
  for (int i = 0; i < static_cast<int>(memerna_dps.size()); ++i) {
    if (memerna_computeds_[0].energy != memerna_computeds_[i].energy ||
        memerna_computeds_[0].energy != memerna_ctd_efns[i] ||
        memerna_computeds_[0].energy != memerna_optimal_efns[i])
      errors.push_back(sfmt("memerna %d: %d (dp) %d (ctd efn) %d (efn) != mfe %d", i,
          memerna_computeds_[i].energy, memerna_ctd_efns[i], memerna_optimal_efns[i],
          memerna_computeds_[0].energy));
  }

  return errors;
}

Error Fuzzer::RnastructureComputeAndCheckState() {
  Error errors;
#ifdef USE_RNASTRUCTURE
  auto rnastructure_computed = rnastructure_->FoldAndDpTable(r_, &rnastructure_dp_);
  auto rnastructure_efn = rnastructure_->Efn(rnastructure_computed.s);
  if (memerna_computeds_[0].energy != rnastructure_computed.energy ||
      memerna_computeds_[0].energy != rnastructure_efn)
    errors.push_back(sfmt("mfe: rnastructure %d (dp), %d (efn) != mfe %d",
        rnastructure_computed.energy, rnastructure_efn, memerna_computeds_[0].energy));
#endif  // USE_RNASTRUCTURE
  return errors;
}

Error Fuzzer::CheckBruteForce() {
  Error errors;
  ModelCfg cfg(ModelCfg::TableAlg::TWO, ModelCfg::SuboptimalAlg::ONE, ModelCfg::PartitionAlg::ZERO);
  Context ctx(r_, em_, cfg);

  if (cfg_.subopt) {
    auto brute_subopt = SuboptimalBruteForce(r_, em_, cfg_.brute_subopt_max);
    auto memerna_subopt = ctx.SuboptimalIntoVector(true, -1, cfg_.brute_subopt_max);

    AppendErrors(
        errors, MaybePrepend(CheckSuboptimalResult(brute_subopt, true), "brute suboptimal:"));
    AppendErrors(
        errors, MaybePrepend(CheckSuboptimalResult(memerna_subopt, true), "memerna suboptimal:"));
    AppendErrors(errors,
        MaybePrepend(CheckSuboptimalResultPair(brute_subopt, memerna_subopt),
            "brute vs memerna suboptimal:"));
  }

  if (cfg_.partition) {
    auto memerna_partition = ctx.Partition();

    auto brute_partition = PartitionBruteForce(r_, em_);
    // Types for the partition function are meant to be a bit configurable, so use sstream here.
    if (!equ(brute_partition.first.q, memerna_partition.q)) {
      std::stringstream sstream;
      sstream << "q: brute partition " << brute_partition.first.q << " != memerna "
              << memerna_partition.q
              << "; difference: " << brute_partition.first.q - memerna_partition.q;
      errors.push_back(sstream.str());
    }

    for (int st = 0; st < static_cast<int>(r_.size()); ++st) {
      for (int en = 0; en < static_cast<int>(r_.size()); ++en) {
        if (!equ(brute_partition.first.p[st][en][0], memerna_partition.p[st][en][0])) {
          std::stringstream sstream;
          sstream << "memerna " << st << " " << en << ": " << memerna_partition.p[st][en][0]
                  << " != brute force " << brute_partition.first.p[st][en][0] << "; difference: "
                  << brute_partition.first.p[st][en][0] - memerna_partition.p[st][en][0];
          errors.push_back(sstream.str());
        }
      }
    }
  }
  return errors;
}

Error Fuzzer::CheckPartition() {
  Error errors;
  std::vector<partition::Partition> memerna_partitions;
  for (auto partition_alg : ModelCfg::PARTITION_ALGS) {
    Context ctx(
        r_, em_, ModelCfg(ModelCfg::TableAlg::TWO, ModelCfg::SuboptimalAlg::ONE, partition_alg));
    memerna_partitions.emplace_back(ctx.Partition());
  }

  for (int i = 0; i < static_cast<int>(memerna_partitions.size()); ++i) {
    if (!equ(memerna_partitions[i].q, memerna_partitions[0].q)) {
      std::stringstream sstream;
      sstream << "q: memerna partition " << i << ": " << memerna_partitions[i].q
              << " != " << memerna_partitions[0].q
              << "; difference: " << memerna_partitions[i].q - memerna_partitions[0].q;
      errors.push_back(sstream.str());
    }

    for (int st = 0; st < static_cast<int>(r_.size()); ++st) {
      for (int en = 0; en < static_cast<int>(r_.size()); ++en) {
        if (!equ(memerna_partitions[i].p[st][en][0], memerna_partitions[0].p[st][en][0])) {
          std::stringstream sstream;
          sstream << "memerna " << i << " at " << st << " " << en << ": "
                  << memerna_partitions[i].p[st][en][0]
                  << " != " << memerna_partitions[0].p[st][en][0] << "; difference: "
                  << memerna_partitions[i].p[st][en][0] - memerna_partitions[0].p[st][en][0];
          errors.push_back(sstream.str());
        }
      }
    }
  }

#ifdef USE_RNASTRUCTURE
  if (cfg_.partition_rnastructure) {
    auto rnastructure_part = rnastructure_->Partition(r_);
    // Types for the partition function are meant to be a bit configurable, so use sstream here.
    if (!equ(rnastructure_part.first.q, memerna_partitions[0].q)) {
      std::stringstream sstream;
      sstream << "q: rnastructure partition " << rnastructure_part.first.q << " != memerna "
              << memerna_partitions[0].q
              << "; difference: " << rnastructure_part.first.q - memerna_partitions[0].q;
      errors.push_back(sstream.str());
    }

    for (int st = 0; st < static_cast<int>(r_.size()); ++st) {
      for (int en = 0; en < static_cast<int>(r_.size()); ++en) {
        if (!equ(rnastructure_part.first.p[st][en][0], memerna_partitions[0].p[st][en][0])) {
          std::stringstream sstream;
          sstream << "memerna " << st << " " << en << ": " << memerna_partitions[0].p[st][en][0]
                  << " != rnastructure " << rnastructure_part.first.p[st][en][0] << "; difference: "
                  << rnastructure_part.first.p[st][en][0] - memerna_partitions[0].p[st][en][0];
          errors.push_back(sstream.str());
        }
      }
    }
  }
#endif  // USE_RNASTRUCTURE

  return errors;
}

}  // namespace mrna::fuzz
