// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "api/bridge/bridge.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/mfe.h"
#include "api/options.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "model/pfn.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::bridge::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_VERBOSE);
  args.RegisterOpt(mrna::OPT_EFN);
  args.RegisterOpt(mrna::OPT_FOLD);
  args.RegisterOpt(mrna::OPT_SUBOPT);
  args.RegisterOpt(mrna::OPT_PFN);
  args.ParseOrExit(argc, argv);

  const bool efn = args.GetOr(mrna::OPT_EFN);
  const bool fold = args.GetOr(mrna::OPT_FOLD);
  const bool subopt = args.GetOr(mrna::OPT_SUBOPT);
  const bool pfn = args.GetOr(mrna::OPT_PFN);

  verify(efn + fold + subopt + pfn == 1, "require exactly one program flag\n{}", args.Usage());

  auto package = mrna::bridge::RnaPackage::FromArgParse(args);
  auto energy_cfg = mrna::erg::EnergyCfg::FromArgParse(args);
  auto pf_cfg = mrna::erg::PseudofreeCfg::FromArgParse(args);
  auto trace_cfg = mrna::trace::TraceCfg::FromArgParse(args);
  auto mfe_alg = args.MaybeGet<mrna::MfeAlg>(mrna::OPT_MFE_ALG);
  auto subopt_alg = args.MaybeGet<mrna::SuboptAlg>(mrna::OPT_SUBOPT_ALG);
  auto pfn_alg = args.MaybeGet<mrna::PfnAlg>(mrna::OPT_PFN_ALG);
  auto subopt_cfg = mrna::subopt::SuboptCfg::FromArgParse(args);

  std::deque<std::string> q(args.Pos().begin(), args.Pos().end());
  const bool read_stdin = q.empty();
  if (efn) {
    while (true) {
      std::string seq;
      std::string db;
      if (read_stdin) {
        getline(std::cin, seq);
        getline(std::cin, db);
        if (!std::cin) break;
      } else {
        if (q.empty()) break;
        verify(q.size() >= 2, "need even number of positional args");
        seq = q.front();
        q.pop_front();
        db = q.front();
        q.pop_front();
      }
      auto [r, s] = mrna::ParseSeqDb(seq, db);
      pf_cfg.Verify(r);
      const auto res = package->Efn(r, s, energy_cfg, pf_cfg, /*given_ctd=*/nullptr,
          /*build_structure=*/args.GetOr(mrna::OPT_VERBOSE));
      fmt::print("{}\n", res.energy);
      if (res.struc)
        for (const auto& line : res.struc->Description()) fmt::print("{}\n", line);
    }
  } else {
    while (true) {
      std::string seq;
      if (read_stdin) {
        getline(std::cin, seq);
        if (!std::cin) break;
      } else {
        if (q.empty()) break;
        seq = q.front();
        q.pop_front();
      }
      auto r = mrna::Primary::FromSeq(seq);
      pf_cfg.Verify(r);

      if (subopt) {
        int64_t strucs = package->Subopt(
            r, mfe_alg, subopt_alg, energy_cfg, pf_cfg,
            [](const mrna::subopt::SuboptResult& c) {
              fmt::print("{} {}\n", c.energy, c.tb.s.ToDb());
            },
            subopt_cfg);
        fmt::print("{} suboptimal structures\n", strucs);
      } else if (fold) {
        const auto res = package->Fold(r, mfe_alg, energy_cfg, pf_cfg, trace_cfg);
        fmt::print("{}\n{}\n", res.mfe.energy, res.tb.s.ToDb());
      } else if (pfn) {
        auto res = package->Pfn(r, pfn_alg, energy_cfg, pf_cfg);
        fmt::print("q: {}\np:\n", res.pfn.q);
        PrintPfn(res.pfn.p);
        fmt::print("\nprobabilities:\n");
        PrintBoltzProbs(res.pfn.prob);
      }
    }
  }
}
