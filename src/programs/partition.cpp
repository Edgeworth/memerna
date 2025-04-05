// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/pfn.h"
#include "model/pfn.h"
#include "model/primary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/util.h"

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

  auto ctx = mrna::Ctx::FromArgParse(args);
  auto res = ctx.Pfn(r);
  fmt::print("q: " FLTFMT "\np:\n", res.pfn.q);
  mrna::PrintPfn(res.pfn.p);
  fmt::print("\nprobabilities:\n");
  mrna::PrintBoltzProbs(res.pfn.prob);

  if (args.Has(mrna::OPT_PFN_PRINT_EXTRA_TABLES)) {
    auto vis = mrna::overloaded{[&](const auto& model) {
      fmt::print("\nHelix probabilities\n");
      mrna::PrintHelixProbs(r, res.pfn, *model);
      fmt::print("\nInner stacking probabilities\n");
      mrna::PrintInnerStackProbs(r, res.pfn, *model);
    }};
    std::visit(vis, ctx.m());
  }
}
