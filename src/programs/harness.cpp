// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#include <iostream>
#include <memory>
#include "bridge/bridge.h"
#include "energy/load_model.h"
#include "parsing.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse({
      {"v", {"be verbose (if possible)"}}, {"e", {"run efn"}}, {"f", {"run fold"}},
      {"subopt-delta", opt_t("maximum energy delta from minimum").Arg("-1")},
  });
  argparse.AddOptions(bridge::BRIDGE_OPTIONS);
  argparse.AddOptions(CONTEXT_OPTIONS);
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  verify_expr(argparse.HasFlag("e") + argparse.HasFlag("f") == 1,
      "require exactly one program flag\n%s", argparse.Usage().c_str());

  const auto package = bridge::RnaPackageFromArgParse(argparse);
  const auto& pos = argparse.GetPositional();
  const bool read_stdin = pos.empty();
  std::deque<std::string> rnaqueue(pos.begin(), pos.end());
  if (argparse.HasFlag("e")) {
    while (1) {
      std::string seq, db;
      if (read_stdin) {
        getline(std::cin, seq);
        getline(std::cin, db);
        if (!std::cin) break;
      } else {
        if (rnaqueue.empty()) break;
        verify_expr(rnaqueue.size() >= 2, "need even number of positional args");
        seq = rnaqueue.front();
        rnaqueue.pop_front();
        db = rnaqueue.front();
        rnaqueue.pop_front();
      }
      const auto secondary = parsing::ParseDotBracketSecondary(seq, db);
      std::string desc;
      const auto res = package->Efn(secondary, argparse.HasFlag("v") ? &desc : nullptr);
      printf("%d\n%s", res, desc.c_str());
    }
  } else {
    while (1) {
      std::string seq;
      if (read_stdin) {
        getline(std::cin, seq);
        if (!std::cin) break;
      } else {
        if (rnaqueue.empty()) break;
        seq = rnaqueue.front();
        rnaqueue.pop_front();
      }
      const auto r = parsing::StringToPrimary(seq);
      int subopt_delta = atoi(argparse.GetOption("subopt-delta").c_str());
      if (subopt_delta >= 0) {
        int num_structures = package->Suboptimal([](const computed_t& c) {
          printf("%d %s\n", c.energy, parsing::PairsToDotBracket(c.s.p).c_str());
        }, r, subopt_delta);
        printf("%d suboptimal structures:\n", num_structures);
      } else {
        const auto res = package->Fold(r);
        printf("%d\n%s\n", res.energy, parsing::PairsToDotBracket(res.s.p).c_str());
      }
    }
  }
}
