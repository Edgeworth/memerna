#include "argparse.h"

namespace memerna {

std::string ArgParse::Parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    bool is_flag = arg[0] == '-';
    while (*arg == '-') ++arg;

    if (is_flag) {
      if (possible_args.count(arg) == 0)
        return sfmt("unknown argument %s", argv[i]);
      if (i + 1 == argc)
        return sfmt("missing argument for flag %s", arg);
      flags[arg] = argv[++i];
    } else {
      positional.push_back(std::move(arg));
    }
  }
  for (const auto& arg : possible_args) {
    if (arg.second.has_arg && flags.count(arg.first) == 0 && !arg.second.has_default)
      return sfmt("missing argument for flag %s", arg.first.c_str());
  }
  return "";
}

std::string ArgParse::Usage() {
  std::string usage = "Flags: \n";
  for (const auto& arg : possible_args) {
    usage += sfmt("  -%s: %s\n", arg.first.c_str(), arg.second.desc.c_str());
  }
  return usage;
}

}
