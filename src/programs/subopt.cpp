// Copyright 2016 E.
#include <cstdio>

#include "compute/energy/load_model.h"
#include "model/context.h"
#include "model/parsing.h"

using mrna::computed_t;
using mrna::context_opt_t;
using mrna::energy_t;
using mrna::opt_t;

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::COMPUTE_ENERGY_OPTIONS);
  args.AddOptions(mrna::CONTEXT_OPTIONS);
  args.AddOptions({{"delta", opt_t("maximum energy delta from minimum").Arg("-1")},
      {"num", opt_t("maximum number of reported structures").Arg("-1")}, {"q", opt_t("quiet")},
      {"sorted", opt_t("if the structures should be sorted")},
      {"ctd-output", opt_t("if we should output CTD data")}});
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  auto opt = ContextOptionsFromArgParse(args);
  opt.table_alg = context_opt_t::TableAlg::TWO;
  mrna::Context ctx(
      mrna::StringToPrimary(pos.front()), mrna::energy::LoadEnergyModelFromArgParse(args), opt);

  const energy_t subopt_delta = atoi(args.GetOption("delta").c_str());
  const int subopt_num = atoi(args.GetOption("num").c_str());
  const bool should_print = !args.HasFlag("q");
  const bool sorted = args.HasFlag("sorted");
  const bool ctd_data = args.HasFlag("ctd-output");
  verify(subopt_delta >= 0 || subopt_num > 0, "nothing to do");

  mrna::subopt::SuboptimalCallback fn = [](const computed_t&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [](const computed_t& c) {
        printf("%d ", c.energy);
        puts(mrna::ComputedToCtdString(c).c_str());
      };
    } else {
      fn = [](const computed_t& c) {
        printf("%d ", c.energy);
        puts(mrna::mfe::internal::grep.c_str());  // meme
      };
    }
  }
  int num_structures = ctx.Suboptimal(fn, sorted, subopt_delta, subopt_num);
  printf("%d suboptimal structures\n", num_structures);
}
