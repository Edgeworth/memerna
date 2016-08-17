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

      if (possible_args[arg].has_arg) {
        if (i + 1 == argc)
          return sfmt("missing argument for flag %s", arg);
        flags[arg] = argv[++i];
      } else {
        flags[arg] = arg;
      }
    } else {
      positional.push_back(std::move(arg));
    }
  }
  for (const auto& argpair : possible_args) {
    const auto& flag = argpair.first;
    const auto& arg = argpair.second;
    verify_expr(!arg.has_arg || arg.has_default, "bad option somehow");
    if (arg.has_arg && flags.count(flag) == 0 && !arg.has_default && arg.required)
      return sfmt("missing argument for flag %s", flag.c_str());
    if (!arg.choices.empty() && arg.choices.count(GetOption(flag)) == 0)
      return sfmt("unrecognised argument for flag %s", flag.c_str());
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
