#ifndef MEMERNA_ARGPARSE_H
#define MEMERNA_ARGPARSE_H

#include "common.h"
#include <unordered_set>
#include <map>

namespace memerna {

class ArgParse {
public:
  ArgParse(const std::map<std::string, std::string>& _possible_args) : possible_args(_possible_args) {}
  ArgParse(const ArgParse&) = delete;
  ArgParse& operator=(const ArgParse&) = delete;

  std::string Parse(int argc, char* argv[]);
  std::string Usage();
  const std::vector<std::string>& GetPositional() const { return positional; }
  bool HasFlag(const std::string& flag) const { return flags.count(flag); }

private:
  std::map<std::string, std::string> possible_args;
  std::unordered_set<std::string> flags;
  std::vector<std::string> positional;
};

}

#endif //MEMERNA_ARGPARSE_H
