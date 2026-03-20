// Copyright 2022 Eliot Courtney.
#include "fuzz/fuzz_cfg.h"
#include "util/argparse.h"
#include "util/error.h"

#ifdef __AFL_FUZZ_TESTCASE_LEN
#include <fmt/core.h>
#include <unistd.h>

#include <boost/json/src.hpp>
#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "fuzz/fuzz_harness.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/enum.h"

__AFL_FUZZ_INIT();

namespace json = boost::json;
namespace {

template <typename T>
std::optional<T> ParseEnumField(const json::object& obj, const char* key) {
  const auto* v = obj.if_contains(key);
  if (!v || !v->is_string()) return std::nullopt;
  auto s = mrna::NormalizeEnumName(std::string(v->as_string()));
  for (const auto& [val, name] : mrna::EnumItems<T>()) {
    if (s == name) return val;
  }
  return std::nullopt;
}

std::optional<uint_fast32_t> ParseSeed(const json::object& obj) {
  const auto* v = obj.if_contains("seed");
  if (!v) return std::nullopt;
  auto r = v->try_to_number<uint_fast32_t>();
  if (r) return *r;
  return std::nullopt;
}

std::optional<std::string> ParseSeq(const json::object& obj) {
  const auto* v = obj.if_contains("seq");
  if (!v || !v->is_string()) return std::nullopt;
  return std::string(v->as_string());
}

mrna::erg::PseudofreeCfg MakePseudofree(const mrna::fuzz::FuzzCfg& fuzz_cfg, std::size_t len) {
  if (!fuzz_cfg.random_pseudofree) return {};

  std::mt19937 pf_rng(fuzz_cfg.seed.value_or(0));
  return {mrna::RandomEnergies(len, mrna::E(-10.0), mrna::E(10.0), pf_rng),
      mrna::RandomEnergies(len, mrna::E(-10.0), mrna::E(10.0), pf_rng)};
}

}  // namespace
#endif

// Limit to 200 by default to avoid OOM.
inline const auto OPT_MAX_LEN = mrna::Opt(mrna::Opt::ARG)
                                    .LongName("max-len")
                                    .Help("limit max length of sequences fuzzed")
                                    .Default(200);

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(OPT_MAX_LEN);
  args.ParseOrExit(argc, argv);

  [[maybe_unused]] const auto max_len = args.Get<int>(OPT_MAX_LEN);

  auto base_cfg = mrna::fuzz::FuzzCfg::FromArgParse(args);

#ifdef __AFL_FUZZ_TESTCASE_LEN
  __AFL_INIT();
  // This must be after __AFL_INIT and before __AFL_LOOP.
  auto buf = reinterpret_cast<const char*>(__AFL_FUZZ_TESTCASE_BUF);
  while (__AFL_LOOP(1000)) {
    const auto len = static_cast<std::size_t>(__AFL_FUZZ_TESTCASE_LEN);

    // Parse JSON testcase.
    boost::system::error_code ec;
    auto jv = json::parse(json::string_view(buf, len), ec);
    if (ec || !jv.is_object()) continue;
    const auto& obj = jv.as_object();

    auto rs = ParseSeq(obj);
    if (!rs || rs->empty()) continue;

    // Copy base config and apply JSON overrides.
    auto fuzz_cfg = base_cfg;
    if (auto seed = ParseSeed(obj)) fuzz_cfg.seed = *seed;
    if (auto em = ParseEnumField<mrna::erg::EnergyModelKind>(obj, "energy_model"))
      fuzz_cfg.energy_model = *em;
    if (auto ctd = ParseEnumField<mrna::erg::EnergyCfg::Ctd>(obj, "ctd"))
      fuzz_cfg.energy_cfg.ctd = *ctd;

    // Parse sequence.
    mrna::Primary seq;
    try {
      if (max_len > 0 && static_cast<int>(rs->size()) > max_len) rs->resize(max_len);
      seq = mrna::Primary::FromSeq(*rs);
    } catch (const std::exception& e) {
      continue;
    }

    auto harness = mrna::fuzz::FuzzHarness(fuzz_cfg, /*should_log=*/false);
    auto pf = MakePseudofree(fuzz_cfg, seq.size());
    const auto res = harness.Run(seq, pf);
    if (!res.empty()) {
      // Print info needed to reproduce with the standalone fuzz program.
      fmt::print("Sequence: {}\n", seq.ToSeq());
      fmt::print("Energy model: {}\n", fuzz_cfg.energy_model);
      fmt::print("Energy cfg: {}\n", fuzz_cfg.energy_cfg);
      if (fuzz_cfg.seed.has_value()) fmt::print("Seed: {}\n", *fuzz_cfg.seed);
      if (!pf.paired.empty()) {
        fmt::print("Pseudofree paired:");
        for (const auto& e : pf.paired) fmt::print(" {}", e);
        fmt::print("\n");
      }
      if (!pf.unpaired.empty()) {
        fmt::print("Pseudofree unpaired:");
        for (const auto& e : pf.unpaired) fmt::print(" {}", e);
        fmt::print("\n");
      }
      for (const auto& s : res) fmt::print("{}\n", s);
      fmt::print("\n");
      abort();
    }
  }
#endif
}
