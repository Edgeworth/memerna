#include "argparse.h"

namespace memerna {

std::string ArgParse::Parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg[0] == '-') {
      auto flag = argv[i] + 1;
      if (possible_args.count(flag) == 0)
        return sfmt("unknown argument %s", argv[i]);
      flags.insert(std::move(flag));
    } else {
      positional.push_back(std::move(arg));
    }
  }
  return "";
}

std::string ArgParse::Usage() {
  std::string usage = "Flags: \n";
  for (const auto& arg : possible_args) {
    usage += sfmt("  -%s: %s\n", arg.first.c_str(), arg.second.c_str());
  }
  return usage;
}

}
