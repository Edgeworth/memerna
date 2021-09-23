// Copyright 2016 Eliot Courtney.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <random>
#include <set>
#include <sstream>

#include "bridge/bridge.h"
#include "bridge/memerna.h"
#include "bridge/rnastructure.h"
#include "common.h"
#include "energy/load_model.h"
#include "energy/structure.h"
#include "fold/brute_fold.h"
#include "fold/fold_constants.h"
#include "fold/fold_globals.h"
#include "parsing.h"

using mrna::opt_t;

namespace mrna {

using fold::PartitionBruteForce;
using fold::SuboptimalBruteForce;
using fold::internal::DP_P;
using fold::internal::DP_SIZE;
using fold::internal::DP_U;
using fold::internal::gdp;

struct cfg_t {
  bool random_model = false;
  bool rnastructure = false;
  uint_fast32_t seed = 0;

  bool subopt = true;
  bool subopt_rnastructure = false;
  int subopt_max = 5000;
  int subopt_delta = 6;  // Same as RNAstructure default.

  int brute_cutoff = 30;
  int brute_subopt_max = 100000;

  bool partition = true;
  bool partition_rnastructure = false;

  std::string Describe() {
    std::string desc;
    if (random_model)
      desc += "random energy models";
    else
      desc += "specified energy model";
    if (subopt) {
      desc += " - testing suboptimal";
      if (subopt_rnastructure) desc += " (including rnastructure)";
      desc += sfmt(" max: %d delta: %d", subopt_max, subopt_delta);
    }
    desc += sfmt(" - brute cutoff: %d subopt-max: %d", brute_cutoff, brute_subopt_max);
    if (partition) {
      desc += " - testing partition";
      if (partition_rnastructure) desc += " (including rnastructure)";
    }
    return desc;
  }
};

const std::map<std::string, opt_t> CFG_OPTIONS = {
    {"random", opt_t("use random energy models (disables comparison to RNAstructure)")},
    {"no-rnastructure", opt_t("disable rnastructure testing")},
    {"no-subopt", opt_t("whether to test suboptimal folding")},
    {"subopt-rnastructure", opt_t("test rnastructure suboptimal folding")},
    {"subopt-max", opt_t("maximum number of substructures for subopt max-delta fuzz").Arg()},
    {"subopt-delta", opt_t("delta for subopt delta fuzz").Arg()},
    {"brute-cutoff", opt_t("maximum rna size to run brute force on").Arg()},
    {"brute-subopt-max", opt_t("maximum number of substructures for brute force fuzz").Arg()},
    {"no-partition", opt_t("whether to test partition function")},
    {"partition-rnastructure", opt_t("test rnastructure partition function")},
};

cfg_t CfgFromArgParse(const ArgParse& args) {
  cfg_t cfg;
  cfg.random_model = args.HasFlag("random");
  cfg.rnastructure = !args.HasFlag("no-rnastructure");
  cfg.subopt = !args.HasFlag("no-subopt");
  cfg.subopt_rnastructure = args.HasFlag("subopt-rnastructure");
  cfg.partition_rnastructure = args.HasFlag("partition-rnastructure");
  if (args.HasFlag("subopt-max")) cfg.subopt_max = atoi(args.GetOption("subopt-max").c_str());
  if (args.HasFlag("subopt-delta")) cfg.subopt_delta = atoi(args.GetOption("subopt-delta").c_str());
  if (args.HasFlag("brute-cutoff")) cfg.brute_cutoff = atoi(args.GetOption("brute-cutoff").c_str());
  if (args.HasFlag("brute-subopt-max"))
    cfg.brute_subopt_max = atoi(args.GetOption("brute-subopt-max").c_str());
  cfg.partition = !args.HasFlag("no-partition");

  verify(!cfg.subopt_rnastructure || cfg.subopt,
      "suboptimal folding testing must be enabled to test rnastructure suboptimal folding");
  verify(!(cfg.random_model && cfg.subopt_rnastructure),
      "cannot use a random energy model with rnastructure");
  verify(cfg.rnastructure || (!cfg.subopt_rnastructure && !cfg.partition_rnastructure),
      "rnastructure must be enabled to use it for suboptimal or partition");
  verify(!cfg.random_model || !cfg.rnastructure,
      "rnastructure testing does not support random models");
  return cfg;
}

inline bool equ(penergy_t a, penergy_t b) { return fabs(a - b) < EP; }

class Fuzzer {
 public:
  typedef std::deque<std::string> error_t;

  Fuzzer(primary_t r_, const cfg_t& cfg_, const energy::EnergyModelPtr em_,
      const bridge::RNAstructure& rnastructure_)
      : N(static_cast<int>(r_.size())), r(std::move(r_)), cfg(cfg_),
        em(cfg.random_model ? energy::LoadRandomEnergyModel(cfg.seed) : em_),
        rnastructure(rnastructure_) {}

  error_t Run() {
    error_t errors;
    AppendErrors(errors, MaybePrependHeader(MemernaComputeAndCheckState(), "memerna:"));
    if (cfg.rnastructure)
      AppendErrors(errors, MaybePrependHeader(RnastructureComputeAndCheckState(), "rnastructure:"));
    AppendErrors(errors, MaybePrependHeader(CheckDpTables(), "dp tables:"));
    if (cfg.subopt) AppendErrors(errors, MaybePrependHeader(CheckSuboptimal(), "suboptimal:"));

    if (static_cast<int>(r.size()) <= cfg.brute_cutoff)
      AppendErrors(errors, MaybePrependHeader(CheckBruteForce(), "brute force:"));

    if (cfg.partition) AppendErrors(errors, MaybePrependHeader(CheckPartition(), "partition:"));

    if (!errors.empty()) {
      if (cfg.random_model)
        errors.push_front(sfmt("Used random energy model with seed: %" PRIuFAST32 "\n", cfg.seed));
      else
        errors.push_front(sfmt("Used specified energy model"));
      errors = MaybePrependHeader(errors,
          sfmt("Difference on len %zu RNA %s:", r.size(), parsing::PrimaryToString(r).c_str()));
    }

    return errors;
  }

 private:
  const int N;
  const primary_t r;
  const cfg_t cfg;
  const energy::EnergyModelPtr em;
  const bridge::RNAstructure& rnastructure;

  // Fuzz state.
  std::vector<computed_t> memerna_computeds;
  std::vector<array3d_t<energy_t, DP_SIZE>> memerna_dps;
  dp_state_t rnastructure_dp;

  error_t MaybePrependHeader(const error_t& main, const std::string& header) {
    if (main.empty()) return main;
    error_t nmain;
    nmain.push_front(header);
    for (auto& error : main) nmain.push_back("  " + error);  // mfw this inefficiency
    return nmain;
  }

  void AppendErrors(error_t& main, error_t&& extra) {
    for (auto& s : extra) main.push_back(std::move(s));
  }

  bool HasDuplicates(const std::vector<computed_t>& computeds) {
    // If energies are different but everything else is the same, it is still a bug.
    std::set<std::pair<secondary_t, std::vector<Ctd>>> suboptimal_set;
    for (const auto& computed : computeds) {
      auto val = std::make_pair(computed.s, computed.base_ctds);
      if (suboptimal_set.count(val)) return true;
      suboptimal_set.insert(val);
    }
    return false;
  }

  error_t CheckSuboptimalResult(const std::vector<computed_t>& subopt, bool has_ctds) {
    error_t errors;
    // Check at least one suboptimal structure.
    if (subopt.empty()) errors.push_back("no structures returned");
    // Check MFE.
    if (!subopt.empty() && memerna_computeds[0].energy != subopt[0].energy)
      errors.push_back(sfmt(
          "lowest structure energy %d != mfe %d", subopt[0].energy, memerna_computeds[0].energy));

    // Only ones with CTDs set can do these tests.
    if (has_ctds) {
      // Check for duplicate structures.
      if (HasDuplicates(subopt)) errors.push_back("has duplicates");

      for (int i = 0; i < static_cast<int>(subopt.size()); ++i) {
        const auto& structure = subopt[i];
        auto suboptimal_efn = energy::ComputeEnergyWithCtds(structure, *em);
        if (suboptimal_efn.energy != structure.energy) {
          errors.push_back(sfmt(
              "structure %d: energy %d != efn %d", i, structure.energy, suboptimal_efn.energy));
          break;
        }

        // Incidentally test ctd parsing.
        auto parsed_computed = parsing::ParseCtdComputed(
            parsing::PrimaryToString(structure.s.r), parsing::ComputedToCtdString(structure));
        parsed_computed.energy = structure.energy;
        if (parsed_computed != structure) {
          errors.push_back(sfmt("structure %d: bug in parsing code", i));
          break;
        }
      }
    }
    return errors;
  }

  error_t CheckSuboptimalResultPair(
      const std::vector<computed_t>& a, const std::vector<computed_t>& b) {
    error_t errors;
    if (a.size() != b.size()) {
      errors.push_back(
          sfmt("first has %zu structures != second has %zu structures", a.size(), b.size()));
    } else {
      for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        if (a[i].energy != b[i].energy) {
          errors.push_back(
              sfmt("structure %d: first %d != second %d", i, a[i].energy, b[i].energy));
          break;
        }
      }
    }
    return errors;
  }

  error_t CheckSuboptimal() {
    error_t errors;
    std::vector<std::vector<computed_t>> memerna_subopts_delta, memerna_subopts_num;
    for (auto subopt_alg : context_opt_t::SUBOPTIMAL_ALGS) {
      context_opt_t options(context_opt_t::TableAlg::TWO, subopt_alg);
      Context ctx(r, em, options);
      memerna_subopts_delta.push_back(ctx.SuboptimalIntoVector(true, cfg.subopt_delta, -1));
      memerna_subopts_num.push_back(ctx.SuboptimalIntoVector(true, -1, cfg.subopt_max));
    }

    for (int i = 0; i < static_cast<int>(memerna_subopts_delta.size()); ++i) {
      AppendErrors(errors,
          MaybePrependHeader(CheckSuboptimalResult(memerna_subopts_delta[i], true),
              sfmt("memerna delta suboptimal %d:", i)));
      AppendErrors(errors,
          MaybePrependHeader(
              CheckSuboptimalResultPair(memerna_subopts_delta[0], memerna_subopts_delta[i]),
              sfmt("memerna 0 vs memerna %d delta suboptimal:", i)));
    }

    for (int i = 0; i < static_cast<int>(memerna_subopts_num.size()); ++i) {
      AppendErrors(errors,
          MaybePrependHeader(CheckSuboptimalResult(memerna_subopts_num[i], true),
              sfmt("memerna num suboptimal %d:", i)));
      AppendErrors(errors,
          MaybePrependHeader(
              CheckSuboptimalResultPair(memerna_subopts_num[0], memerna_subopts_num[i]),
              sfmt("memerna 0 vs memerna %d num suboptimal:", i)));
    }

    if (cfg.subopt_rnastructure) {
      // Suboptimal folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does
      // strange things
      // when the energy for suboptimal structures is 0 or above.
      if (memerna_computeds[0].energy < -cfg.subopt_delta) {
        const auto rnastructure_subopt = rnastructure.SuboptimalIntoVector(r, cfg.subopt_delta);
        AppendErrors(errors,
            MaybePrependHeader(
                CheckSuboptimalResult(rnastructure_subopt, false), "rnastructure suboptimal:"));
        AppendErrors(errors,
            MaybePrependHeader(
                CheckSuboptimalResultPair(memerna_subopts_delta[0], rnastructure_subopt),
                "memerna vs rnastructure suboptimal:"));
      }
    }

    return errors;
  }

  error_t CheckDpTables() {
    error_t errors;
    for (int st = N - 1; st >= 0; --st) {
      for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
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
          if (cfg.rnastructure && (a == DP_P || a == DP_U)) {
            energy_t rnastructureval = a == DP_P ? rnastructure_dp.v.f(st + 1, en + 1)
                                                 : rnastructure_dp.w.f(st + 1, en + 1);
            if (((memerna0 < CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
                    (memerna0 < CAP_E && memerna0 != rnastructureval))) {
              errors.push_back(sfmt("rnastructure at %d %d %d: %d != %d", st, en, a,
                  rnastructureval, memerna_dps[0][st][en][a]));
              goto loopend;
            }
          }
        }
      }
    }
  loopend:
    return errors;
  }

  error_t MemernaComputeAndCheckState() {
    error_t errors;
    // Memerna.
    std::vector<energy_t> memerna_ctd_efns;
    std::vector<energy_t> memerna_optimal_efns;
    for (auto table_alg : context_opt_t::TABLE_ALGS) {
      Context ctx(r, em, context_opt_t(table_alg));
      auto computed = ctx.Fold();
      memerna_dps.emplace_back(std::move(gdp));
      // First compute with the CTDs that fold returned to check the energy.
      memerna_ctd_efns.push_back(energy::ComputeEnergyWithCtds(computed, *em).energy);
      // Also check that the optimal CTD configuration has the same energy.
      // Note that it might not be the same, so we can't do an equality check.
      memerna_optimal_efns.push_back(energy::ComputeEnergy(computed.s, *em).energy);
      memerna_computeds.push_back(std::move(computed));
    }

    // Check memerna energies.
    for (int i = 0; i < static_cast<int>(memerna_dps.size()); ++i) {
      if (memerna_computeds[0].energy != memerna_computeds[i].energy ||
          memerna_computeds[0].energy != memerna_ctd_efns[i] ||
          memerna_computeds[0].energy != memerna_optimal_efns[i])
        errors.push_back(sfmt("memerna %d: %d (dp) %d (ctd efn) %d (efn) != mfe %d", i,
            memerna_computeds[i].energy, memerna_ctd_efns[i], memerna_optimal_efns[i],
            memerna_computeds[0].energy));
    }

    return errors;
  }

  error_t RnastructureComputeAndCheckState() {
    error_t errors;
    auto rnastructure_computed = rnastructure.FoldAndDpTable(r, &rnastructure_dp);
    auto rnastructure_efn = rnastructure.Efn(rnastructure_computed.s);
    if (memerna_computeds[0].energy != rnastructure_computed.energy ||
        memerna_computeds[0].energy != rnastructure_efn)
      errors.push_back(sfmt("mfe: rnastructure %d (dp), %d (efn) != mfe %d",
          rnastructure_computed.energy, rnastructure_efn, memerna_computeds[0].energy));
    return errors;
  }

  error_t CheckBruteForce() {
    error_t errors;
    context_opt_t options(context_opt_t::TableAlg::TWO, context_opt_t::SuboptimalAlg::ONE,
        context_opt_t::PartitionAlg::ZERO);
    Context ctx(r, em, options);

    if (cfg.subopt) {
      auto brute_subopt = SuboptimalBruteForce(r, *em, cfg.brute_subopt_max);
      auto memerna_subopt = ctx.SuboptimalIntoVector(true, -1, cfg.brute_subopt_max);

      AppendErrors(errors,
          MaybePrependHeader(CheckSuboptimalResult(brute_subopt, true), "brute suboptimal:"));
      AppendErrors(errors,
          MaybePrependHeader(CheckSuboptimalResult(memerna_subopt, true), "memerna suboptimal:"));
      AppendErrors(errors,
          MaybePrependHeader(CheckSuboptimalResultPair(brute_subopt, memerna_subopt),
              "brute vs memerna suboptimal:"));
    }

    if (cfg.partition) {
      auto memerna_partition = ctx.Partition();

      auto brute_partition = PartitionBruteForce(r, *em);
      // Types for the partition function are meant to be a bit configurable, so use sstream here.
      if (!equ(brute_partition.first.q, memerna_partition.q)) {
        std::stringstream sstream;
        sstream << "q: brute partition " << brute_partition.first.q << " != memerna "
                << memerna_partition.q
                << "; difference: " << brute_partition.first.q - memerna_partition.q;
        errors.push_back(sstream.str());
      }

      for (int st = 0; st < N; ++st) {
        for (int en = 0; en < N; ++en) {
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

  error_t CheckPartition() {
    error_t errors;
    std::vector<partition::partition_t> memerna_partitions;
    for (auto partition_alg : context_opt_t::PARTITION_ALGS) {
      Context ctx(r, em,
          context_opt_t(
              context_opt_t::TableAlg::TWO, context_opt_t::SuboptimalAlg::ONE, partition_alg));
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

      for (int st = 0; st < N; ++st) {
        for (int en = 0; en < N; ++en) {
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

    if (cfg.partition_rnastructure) {
      auto rnastructure_part = rnastructure.Partition(r);
      // Types for the partition function are meant to be a bit configurable, so use sstream here.
      if (!equ(rnastructure_part.first.q, memerna_partitions[0].q)) {
        std::stringstream sstream;
        sstream << "q: rnastructure partition " << rnastructure_part.first.q << " != memerna "
                << memerna_partitions[0].q
                << "; difference: " << rnastructure_part.first.q - memerna_partitions[0].q;
        errors.push_back(sstream.str());
      }

      for (int st = 0; st < N; ++st) {
        for (int en = 0; en < N; ++en) {
          if (!equ(rnastructure_part.first.p[st][en][0], memerna_partitions[0].p[st][en][0])) {
            std::stringstream sstream;
            sstream << "memerna " << st << " " << en << ": " << memerna_partitions[0].p[st][en][0]
                    << " != rnastructure " << rnastructure_part.first.p[st][en][0]
                    << "; difference: "
                    << rnastructure_part.first.p[st][en][0] - memerna_partitions[0].p[st][en][0];
            errors.push_back(sstream.str());
          }
        }
      }
    }
    return errors;
  }
};

}  // namespace mrna

#ifdef __AFL_FUZZ_TESTCASE_LEN
#include <unistd.h>  // For __AFL_FUZZ_TESTCASE_LEN

__AFL_FUZZ_INIT();
#endif

int main(int argc, char* argv[]) {
  std::mt19937 eng(uint_fast32_t(time(nullptr)));
  mrna::ArgParse args({
      {"print-interval", opt_t("status update every n seconds").Arg("-1")},
      {"afl", opt_t("reads one rna from stdin and fuzzes - useful for use with afl")},
  });
  args.AddOptions(mrna::CFG_OPTIONS);
  args.AddOptions(mrna::energy::ENERGY_OPTIONS);
  args.ParseOrExit(argc, argv);

  const mrna::bridge::RNAstructure rnastructure(args.GetOption("data-path"), false);
  const auto em = mrna::energy::LoadEnergyModelFromArgParse(args);

  auto cfg = CfgFromArgParse(args);
  const bool afl_mode = args.HasFlag("afl");
  verify(!cfg.rnastructure || !args.HasFlag("seed"),
      "seed option incompatible with rnastructure testing");

  if (afl_mode) {
#ifdef __AFL_FUZZ_TESTCASE_LEN
    __AFL_INIT();
    while (__AFL_LOOP(1000)) {
      std::string data(
          reinterpret_cast<const char*>(__AFL_FUZZ_TESTCASE_BUF), __AFL_FUZZ_TESTCASE_LEN);
      if (data.size() > 0) {
        cfg.seed = eng();
        mrna::Fuzzer fuzzer(mrna::parsing::StringToPrimary(data), cfg, em, rnastructure);
        const auto res = fuzzer.Run();
        if (!res.empty()) abort();
      }
    }
#endif
  } else {
    auto pos = args.GetPositional();
    verify(pos.size() == 2, "require min and max length");
    const int min_len = atoi(pos[0].c_str());
    const int max_len = atoi(pos[1].c_str());
    const auto interval = atoi(args.GetOption("print-interval").c_str());

    verify(min_len > 0, "invalid min length");
    verify(max_len >= min_len, "invalid max len");
    std::uniform_int_distribution<int> len_dist(min_len, max_len);

    printf("Fuzzing [%d, %d] len RNAs - %s\n", min_len, max_len, cfg.Describe().c_str());

    // Normal mode.
    auto start_time = std::chrono::steady_clock::now();
    for (int64_t i = 0;; ++i) {
      if (interval > 0 &&
          std::chrono::duration_cast<std::chrono::seconds>(
              std::chrono::steady_clock::now() - start_time)
                  .count() > interval) {
        printf("Fuzzed %" PRId64 " RNA\n", i);
        start_time = std::chrono::steady_clock::now();
      }
      int len = len_dist(eng);
      auto r = mrna::GenerateRandomPrimary(len, eng);

      cfg.seed = eng();
      mrna::Fuzzer fuzzer(r, cfg, em, rnastructure);
      const auto res = fuzzer.Run();
      if (!res.empty()) {
        for (const auto& s : res) printf("%s\n", s.c_str());
        printf("\n");
      }
    }
  }
}
