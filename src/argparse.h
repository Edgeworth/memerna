#ifndef MEMERNA_ARGPARSE_H
#define MEMERNA_ARGPARSE_H

#include <map>
#include <unordered_map>
#include <unordered_set>
#include "common.h"

namespace memerna {

class ArgParse {
public:
  struct option_t {
    option_t(const std::string& desc_ = "", const std::string& default_arg_ = "",
        bool has_default_ = false, bool has_arg_ = false, bool required_ = false)
        : desc(desc_), default_arg(default_arg_), choices(), has_default(has_default_),
          has_arg(has_arg_), required(required_) {}

    option_t& Arg(const std::string& default_, const std::unordered_set<std::string>& choices_) {
      default_arg = default_;
      choices = choices_;
      has_default = true;
      has_arg = true;
      return *this;
    }

    option_t& Arg(const std::string& default_) { return Arg(default_, {}); }

    option_t& Arg(const std::unordered_set<std::string>& choices_) {
      choices = choices_;
      has_arg = true;
      return *this;
    }

    option_t& Arg() {
      has_arg = true;
      return *this;
    }

    option_t& Require() {
      required = true;
      return *this;
    }

    std::string Desc() const;

    std::string desc;
    std::string default_arg;
    std::unordered_set<std::string> choices;
    bool has_default;
    bool has_arg;
    bool required;
  };

  ArgParse(const std::map<std::string, option_t>& possible_args_) : possible_args(possible_args_) {}

  ArgParse() = default;
  ArgParse(const ArgParse&) = delete;
  ArgParse& operator=(const ArgParse&) = delete;

  void AddOptions(const std::map<std::string, option_t>& possible_args_);
  std::string Parse(int argc, char* argv[]);
  void ParseOrExit(int argc, char* argv[]);
  std::string Usage() const;

  const std::vector<std::string>& GetPositional() const { return positional; }

  bool HasFlag(const std::string& flag) const {
    if (possible_args.count(flag) && possible_args.find(flag)->second.has_default) return true;
    return flags.count(flag) > 0;
  }

  std::string GetOption(const std::string& flag) const {
    auto flagiter = flags.find(flag);
    auto positer = possible_args.find(flag);
    if (flagiter != flags.end()) return flagiter->second;
    if (positer != possible_args.end()) return positer->second.default_arg;
    return "";
  }

private:
  std::map<std::string, option_t> possible_args;
  std::unordered_map<std::string, std::string> flags;
  std::vector<std::string> positional;
};
}

#endif  // MEMERNA_ARGPARSE_H
