// Copyright 2016 Eliot Courtney.
#ifndef UTIL_ARGPARSE_H_
#define UTIL_ARGPARSE_H_

#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "util/error.h"

namespace mrna {

struct Opt {
  Opt(std::string desc_ = "", std::string default_arg_ = "", bool has_default_ = false,
      bool has_arg_ = false, bool required_ = false)
      : desc(std::move(desc_)), default_arg(std::move(default_arg_)), choices(),
        has_default(has_default_), has_arg(has_arg_), required(required_) {}

  Opt& Arg(std::string default_arg_, std::unordered_set<std::string> choices_) {
    default_arg = std::move(default_arg_);
    choices = std::move(choices_);
    has_default = true;
    has_arg = true;
    return *this;
  }

  Opt& Arg(std::string default_arg_) { return Arg(std::move(default_arg_), {}); }

  Opt& Arg(std::unordered_set<std::string> choices_) {
    choices = std::move(choices_);
    has_arg = true;
    return *this;
  }

  Opt& Arg() {
    has_arg = true;
    return *this;
  }

  Opt& Require() {
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

class ArgParse {
 public:
  explicit ArgParse(std::map<std::string, Opt> possible_args)
      : possible_args_(std::move(possible_args)) {}

  ArgParse() = default;
  ArgParse(const ArgParse&) = delete;
  ArgParse& operator=(const ArgParse&) = delete;

  void AddOptions(const std::map<std::string, Opt>& possible_args_);
  std::string Parse(int argc, char* argv[]);
  void ParseOrExit(int argc, char* argv[]);
  std::string Usage() const;

  const std::vector<std::string>& GetPositional() const { return positional_; }

  bool HasFlag(const std::string& flag) const {
    if (possible_args_.count(flag) && possible_args_.find(flag)->second.has_default) return true;
    return flags.count(flag) > 0;
  }

  std::string GetOption(const std::string& flag) const {
    auto flagiter = flags.find(flag);
    auto positer = possible_args_.find(flag);
    if (flagiter != flags.end()) return flagiter->second;
    if (positer != possible_args_.end()) return positer->second.default_arg;
    error("missing flag %s", flag.c_str());
    return "";
  }

 private:
  std::map<std::string, Opt> possible_args_;
  std::unordered_map<std::string, std::string> flags;
  std::vector<std::string> positional_;
};

}  // namespace mrna

#endif  // UTIL_ARGPARSE_H_
