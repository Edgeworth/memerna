// Copyright 2016 Eliot Courtney.

#include <fmt/base.h>

#include "api/ctx/algorithm.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/pfn.h"
#include "model/pfn.h"
#include "model/primary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/float.h"

namespace mrna {

inline const Opt OPT_PFN_PRINT_EXTRA_TABLES =
    Opt(Opt::FLAG).LongName("print-extra-tables").Help("whether to print extra tables in pfn");

}  // namespace mrna

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_PFN_PRINT_EXTRA_TABLES);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromSeq(args.Pos(0));
  const auto cfg = mrna::erg::EnergyCfgPfn::FromArgParse(args);
  auto pf = mrna::erg::PseudofreeCfg::FromArgParse(args);
  pf.Verify(r);

  auto ctx = mrna::Ctx::FromArgParse(args);
  auto pfn_alg = args.MaybeGet<mrna::PfnAlg>(mrna::OPT_PFN_ALG);
  auto res = ctx.Pfn(r, pfn_alg, cfg, pf);
  fmt::print("q: " FLTFMT "\np:\n", res.pfn.q);
  mrna::PrintPfn(res.pfn.p);
  fmt::print("\nprobabilities:\n");
  mrna::PrintBoltzProbs(res.pfn.prob);

  if (args.GetOr(mrna::OPT_PFN_PRINT_EXTRA_TABLES)) {
    mrna::PrintExtraPfnTables(r, res.pfn, ctx.BackendForPfn(pfn_alg, cfg, pf).m, pf);
  }
}
