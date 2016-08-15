#ifndef MEMERNA_ARGPARSE_H
#define MEMERNA_ARGPARSE_H

#include "common.h"
#include <unordered_set>
#include <map>
#include <unordered_map>

namespace memerna {



// TODO: requirements, etc
class ArgParse {
public:
  struct option_t {
    option_t() = default;
    option_t(const std::string& _desc) : desc(_desc) {}
    option_t& Arg(const std::string& _default) {
      default_arg = _default;
      has_default = true;
      has_arg = true;
      return *this;
    }
    option_t& Arg() {
      has_arg = true;
      return *this;
    }
    std::string desc;
    std::string default_arg;
    bool has_default;
    bool has_arg;
  };

  ArgParse(const std::map<std::string, option_t>& _possible_args) : possible_args(_possible_args) {}
  ArgParse(const ArgParse&) = delete;
  ArgParse& operator=(const ArgParse&) = delete;

  std::string Parse(int argc, char* argv[]);
  std::string Usage();

  const std::vector<std::string>& GetPositional() const {return positional;}

  bool HasFlag(const std::string& flag) const {return bool(flags.count(flag));}
  std::string GetOption(const std::string& flag) {
    if (HasFlag(flag))
      return flags[flag];
    return possible_args[flag].default_arg;
  }

private:
  std::map<std::string, option_t> possible_args;
  std::unordered_map<std::string, std::string> flags;
  std::vector<std::string> positional;
};

}

#endif //MEMERNA_ARGPARSE_H
