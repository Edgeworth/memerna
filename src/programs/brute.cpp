// Copyright 2022 Eliot Courtney.

#include "backends/brute/brute.h"

#include <fmt/core.h>

#include <algorithm>
#include <string>

#include "api/brute/brute_cfg.h"
#include "api/ctx/algorithm.h"
#include "api/ctx/backend.h"
#include "api/ctx/backend_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/options.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "model/secondary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::brute::RegisterOpts(&args);
  mrna::erg::RegisterOptsPseudofree(&args);
  args.ParseOrExit(argc, argv);
  auto backend_cfg = mrna::BackendCfg::FromArgParse(args);
  auto energy_cfg = mrna::erg::EnergyCfg::FromArgParse(args);
  auto pf = mrna::erg::PseudofreeCfg::FromArgParse(args);
  std::string log;
  auto resolved = mrna::ResolveEfn(
      args.MaybeGet<mrna::BackendKind>(mrna::OPT_BACKEND), backend_cfg, energy_cfg, pf, &log);
  if (!resolved.has_value()) fatal("no backend supports this configuration:\n{}", log);
  const auto m = mrna::BackendFromBackendCfg(resolved->backend, backend_cfg);
  auto cfg = mrna::brute::BruteCfg::FromArgParse(args);
  if (cfg.mfe) {
    cfg.subopt = true;
    cfg.subopt_cfg.strucs = std::max(cfg.subopt_cfg.strucs, 1);
    cfg.subopt_cfg.sorted = true;
  }

  verify(args.PosSize() == 1, "requires primary sequence");
  auto r = mrna::Primary::FromSeq(args.Pos(0));
  pf.Verify(r);
  auto res = mrna::md::brute::Brute(r, m, backend_cfg, energy_cfg, pf, cfg).Run();

  if (args.GetOr(mrna::OPT_FOLD)) {
    const auto& mfe = *res.subopts.begin();
    fmt::print("{}\n", mfe.energy);
    fmt::print("{}\n", mfe.tb.s.ToDb());
    fmt::print("{}\n", energy_cfg.ToCtdString(mfe.tb.s, mfe.tb.ctd));
  }

  if (args.GetOr(mrna::OPT_SUBOPT)) {
    for (const auto& s : res.subopts) {
      fmt::print("{}\n", s.energy);
      fmt::print("{}\n", s.tb.s.ToDb());
      fmt::print("{}\n", energy_cfg.ToCtdString(s.tb.s, s.tb.ctd));
    }
  }

  if (args.GetOr(mrna::OPT_PFN)) {
    fmt::print("q: {}\n", res.pfn.q);
    fmt::print("p:\n");
    mrna::PrintPfn(res.pfn.p);
    fmt::print("\nprobabilities:\n");
    mrna::PrintBoltzProbs(res.pfn.prob);
  }
}
