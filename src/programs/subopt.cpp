// Copyright 2016 E.
#include <cstdio>

#include "compute/energy/load_model.h"
#include "model/context.h"

using mrna::Computed;
using mrna::Energy;
using mrna::ModelCfg;
using mrna::Opt;

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::ENERGY_OPTS);
  args.AddOptions(mrna::MODEL_OPTS);
  args.AddOptions({{"delta", Opt("maximum energy delta from minimum").Arg("-1")},
      {"num", Opt("maximum number of reported structures").Arg("-1")}, {"q", Opt("quiet")},
      {"sorted", Opt("if the structures should be sorted")},
      {"ctd-output", Opt("if we should output CTD data")}});
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  auto cfg = ModelCfgFromArgParse(args);
  cfg.table_alg = ModelCfg::TableAlg::TWO;
  mrna::Context ctx(
      mrna::StringToPrimary(pos.front()), mrna::energy::LoadEnergyModelFromArgParse(args), cfg);

  const Energy subopt_delta = atoi(args.GetOption("delta").c_str());
  const int subopt_num = atoi(args.GetOption("num").c_str());
  const bool should_print = !args.HasFlag("q");
  const bool sorted = args.HasFlag("sorted");
  const bool ctd_data = args.HasFlag("ctd-output");
  verify(subopt_delta >= 0 || subopt_num > 0, "nothing to do");

  mrna::subopt::SuboptimalCallback fn = [](const Computed&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [](const Computed& c) {
        printf("%d ", c.energy);
        puts(mrna::ComputedToCtdString(c).c_str());
      };
    } else {
      fn = [](const Computed& c) {
        printf("%d ", c.energy);
        puts(mrna::PairsToDotBracket(c.s.p).c_str());  // meme
      };
    }
  }
  int num_structures = ctx.Suboptimal(fn, sorted, subopt_delta, subopt_num);
  printf("%d suboptimal structures\n", num_structures);
}
