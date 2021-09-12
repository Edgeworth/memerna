// Copyright 2016 E.
#include <cstdio>

#include "context.h"
#include "energy/load_model.h"
#include "parsing.h"

using namespace mrna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(CONTEXT_OPTIONS);
  argparse.AddOptions({{"delta", opt_t("maximum energy delta from minimum").Arg("-1")},
      {"num", opt_t("maximum number of reported structures").Arg("-1")}, {"q", opt_t("quiet")},
      {"sorted", opt_t("if the structures should be sorted")},
      {"ctd-output", opt_t("if we should output CTD data")}});
  argparse.ParseOrExit(argc, argv);
  const auto& pos = argparse.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  auto opt = ContextOptionsFromArgParse(argparse);
  opt.table_alg = context_opt_t::TableAlg::TWO;
  Context ctx(
      parsing::StringToPrimary(pos.front()), energy::LoadEnergyModelFromArgParse(argparse), opt);

  const energy_t subopt_delta = atoi(argparse.GetOption("delta").c_str());
  const int subopt_num = atoi(argparse.GetOption("num").c_str());
  const bool should_print = !argparse.HasFlag("q");
  const bool sorted = argparse.HasFlag("sorted");
  const bool ctd_data = argparse.HasFlag("ctd-output");
  verify(subopt_delta >= 0 || subopt_num > 0, "nothing to do");

  fold::SuboptimalCallback fn = [](const computed_t&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [](const computed_t& c) {
        printf("%d ", c.energy);
        puts(parsing::ComputedToCtdString(c).c_str());
      };
    } else {
      fn = [](const computed_t& c) {
        printf("%d ", c.energy);
        puts(fold::internal::grep.c_str());  // meme
      };
    }
  }
  int num_structures = ctx.Suboptimal(fn, sorted, subopt_delta, subopt_num);
  printf("%d suboptimal structures\n", num_structures);
}
