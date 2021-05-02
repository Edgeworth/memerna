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
#include "argparse.h"

namespace memerna {

std::string opt_t::Desc() const {
  auto res = desc;
  if (has_default) res += sfmt(" [%s]", default_arg.c_str());
  if (!choices.empty()) {
    res += " (";
    for (auto iter = choices.begin(); iter != choices.end(); ++iter) {
      if (iter != choices.begin()) res += ", ";
      res += *iter;
    }
    res += ")";
  }
  return res;
}

std::string ArgParse::Parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    const bool is_flag = arg[0] == '-';
    while (*arg == '-') ++arg;

    if (is_flag) {
      if (possible_args.count(arg) == 0) return sfmt("unknown argument %s", argv[i]);

      if (possible_args[arg].has_arg) {
        if (i + 1 == argc) return sfmt("missing argument for flag %s", arg);
        flags[arg] = argv[++i];
      } else {
        flags[arg] = arg;
      }
    } else {
      positional.push_back(arg);
    }
  }
  for (const auto& argpair : possible_args) {
    const auto& flag = argpair.first;
    const auto& arg = argpair.second;
    verify_expr(!arg.has_default || arg.has_arg, "bad option somehow");
    if (arg.has_arg && flags.count(flag) == 0 && !arg.has_default && arg.required)
      return sfmt("missing argument for flag %s", flag.c_str());
    if (!arg.choices.empty() && arg.choices.count(GetOption(flag)) == 0)
      return sfmt("unrecognised argument for flag %s", flag.c_str());
  }
  return "";
}

std::string ArgParse::Usage() const {
  std::string usage = "Flags: \n";
  for (const auto& arg : possible_args) {
    usage += sfmt("  -%s: %s\n", arg.first.c_str(), arg.second.Desc().c_str());
  }
  return usage;
}

void ArgParse::AddOptions(const std::map<std::string, opt_t>& possible_args_) {
  for (const auto& argpair : possible_args_) {
    verify_expr(
        possible_args.count(argpair.first) == 0, "duplicate argument %s", argpair.first.c_str());
    possible_args.insert(argpair);
  }
}

void ArgParse::ParseOrExit(int argc, char** argv) {
  const auto ret = Parse(argc, argv);
  verify_expr(ret.size() == 0, "%s\n%s\n", ret.c_str(), Usage().c_str());
}
}  // namespace memerna
