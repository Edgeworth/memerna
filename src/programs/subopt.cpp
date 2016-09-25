#include <cstdio>
#include "energy/load_model.h"
#include "fold/context.h"
#include "parsing.h"

using namespace memerna;



int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.AddOptions({
      {"delta", ArgParse::option_t("maximum energy delta from minimum").Arg("-1")},
      {"num", ArgParse::option_t("maximum number of reported structures").Arg("-1")},
      {"q", ArgParse::option_t("quiet")},
      {"sorted", ArgParse::option_t("if the structures should be sorted")}
  });
  argparse.ParseOrExit(argc, argv);
  const auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  auto opt = fold::ContextOptionsFromArgParse(argparse);
  opt.table_alg = fold::context_options_t::TableAlg::TWO;
  fold::Context ctx(
      parsing::StringToPrimary(pos.front()), energy::LoadEnergyModelFromArgParse(argparse), opt);

  const energy_t subopt_delta = atoi(argparse.GetOption("delta").c_str());
  const int subopt_num = atoi(argparse.GetOption("num").c_str());
  const bool should_print = !argparse.HasFlag("q");
  const bool sorted = argparse.HasFlag("sorted");
  verify_expr(subopt_delta >= 0 || subopt_num > 0, "nothing to do");

  int num_structures = 0;
  if (should_print) {
    num_structures = ctx.Suboptimal([](const computed_t& c) {
      printf("%d %s\n", c.energy, parsing::ComputedToCtdString(c).c_str());
    }, sorted, subopt_delta, subopt_num);
  } else {
    num_structures = ctx.Suboptimal([](const computed_t& c) {}, sorted, subopt_delta, subopt_num);
  }
  printf("%d suboptimal structures\n", num_structures);
}
