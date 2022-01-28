// Copyright 2016 Eliot Courtney.
#ifndef UTIL_ARGPARSE_H_
#define UTIL_ARGPARSE_H_

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "util/error.h"

namespace mrna {

struct Opt {
 public:
  Opt& LongName(std::string n) {
    longname_ = std::move(n);
    return *this;
  }

  Opt& ShortName(std::string n) {
    shortname_ = std::move(n);
    return *this;
  }

  Opt& Arg() {
    has_arg_ = true;
    return *this;
  }

  Opt& Choice(std::set<std::string> c) {
    choices_ = std::move(c);
    has_arg_ = true;
    return *this;
  }

  Opt& Default(std::string d) {
    has_arg_ = true;
    has_default_ = true;
    default_arg_ = std::move(d);
    return *this;
  }

  constexpr Opt& Required() {
    required_ = true;
    return *this;
  }

  Opt& Help(std::string h) {
    help_ = std::move(h);
    return *this;
  }

  std::string Desc() const;

  auto operator<=>(const Opt&) const = default;

  const std::string& longname() const { return longname_; }
  const std::string& shortname() const { return shortname_; }
  const std::string& help() const { return help_; }
  const std::string& default_arg() const { return default_arg_; }
  const std::set<std::string>& choices() const { return choices_; }
  bool has_default() const { return has_default_; }
  bool has_arg() const { return has_arg_; }
  bool required() const { return required_; }

 private:
  std::string longname_;
  std::string shortname_;
  std::string help_;
  std::string default_arg_;
  std::set<std::string> choices_;
  bool has_default_;
  bool has_arg_;
  bool required_;
};

// Some common options
inline const Opt OPT_VERBOSE = Opt().LongName("verbose").ShortName("v").Help("verbose output");
inline const Opt OPT_QUIET = Opt().LongName("quiet").ShortName("q").Help("quiet output");
inline const Opt OPT_AFL = mrna::Opt().LongName("afl").Default("-1").Help("run in afl-fuzz mode");

class ArgParse {
 public:
  ArgParse() = default;
  ArgParse(const ArgParse&) = delete;
  ArgParse& operator=(const ArgParse&) = delete;

  void RegisterOpt(const Opt& opt);
  std::string Usage() const;

  std::string Parse(int argc, char* argv[]);
  void ParseOrExit(int argc, char* argv[]);

  const Opt& Lookup(const std::string& name) const;
  bool Has(const Opt& opt) const;
  std::string Get(const Opt& opt) const;

  const std::vector<std::string>& positional() const { return positional_; }

 private:
  std::map<std::string, Opt> longname_;
  std::map<std::string, Opt> shortname_;
  std::set<Opt> opts_;
  std::map<Opt, std::string> values_;
  std::vector<std::string> positional_;
};

}  // namespace mrna

#endif  // UTIL_ARGPARSE_H_
