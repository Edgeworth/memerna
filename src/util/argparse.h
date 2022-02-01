// Copyright 2016 E.
#ifndef UTIL_ARGPARSE_H_
#define UTIL_ARGPARSE_H_

#include <compare>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "util/error.h"
#include "util/string.h"

namespace mrna {

struct Opt {
 public:
  enum Kind {
    FLAG,  // A flag that can be set to true or false, via --flag or --no-flag.
    ARG,  // An argument that can be set to a value.
  };

  Opt() = delete;
  explicit Opt(Kind kind) : kind_(kind) {}

  Opt& LongName(std::string n) {
    longname_ = std::move(n);
    return *this;
  }

  Opt& ShortName(std::string n) {
    shortname_ = std::move(n);
    return *this;
  }

  Opt& Choice(std::set<std::string> c) {
    verify(kind_ == ARG, "cannot set choices for flag");
    choices_ = std::move(c);
    return *this;
  }

  template <typename T>
  Opt& Default(const T& d) {
    has_default_ = true;
    default_ = Conv(d);
    return *this;
  }

  constexpr Opt& Required() {
    required_ = true;
    return *this;
  }

  // Not displayed in command line usage.
  constexpr Opt& Hidden() {
    hidden_ = true;
    return *this;
  }

  Opt& Help(std::string h) {
    help_ = std::move(h);
    return *this;
  }

  std::string Desc() const;

  bool IsInverted() const { return longname_.rfind("no-", 0) == 0; }

  auto operator<=>(const Opt&) const = default;

  Kind kind() const { return kind_; }
  const std::string& longname() const { return longname_; }
  const std::string& shortname() const { return shortname_; }
  const std::string& help() const { return help_; }
  const std::string& default_arg() const { return default_; }
  const std::set<std::string>& choices() const { return choices_; }

  bool has_default() const { return has_default_; }
  bool required() const { return required_; }
  bool hidden() const { return hidden_; }

 private:
  Kind kind_;
  std::string longname_;
  std::string shortname_;
  std::string help_;
  std::string default_;
  std::set<std::string> choices_;
  bool has_default_ = false;
  bool required_ = false;
  bool hidden_ = false;
};

class ArgParse {
 public:
  ArgParse() = default;
  ArgParse(const ArgParse&) = delete;
  ArgParse& operator=(const ArgParse&) = delete;

  ArgParse(ArgParse&&) = default;
  ArgParse& operator=(ArgParse&&) = default;

  void RegisterOpt(const Opt& opt);
  std::string Usage() const;

  std::string Parse(int argc, char* argv[]);
  void ParseOrExit(int argc, char* argv[]);

  const Opt& Lookup(const std::string& name) const;
  const std::vector<std::string>& Pos() const { return pos_; }

  // Note that flags set to false via --no-flag will have Has return true.
  // It returns whether or not this option was specified. Use GetOr instead
  // for flags.
  bool Has(const Opt& opt) const;

  // Prefer this for setting variables over Get<T>(), if those variables
  // have default values.
  template <typename T>
  void MaybeSet(const Opt& opt, T* val) const {
    if (auto iter = values_.find(opt); iter != values_.end())
      *val = Conv<T>(iter->second);  // NOLINT
  }

  // Useful mainly with flags where not specifying them means false (or whatever you pass as the
  // default).
  template <typename T = bool>
  T GetOr(const Opt& opt, T def = T()) const {
    if (auto iter = values_.find(opt); iter != values_.end())
      return Conv<T>(iter->second);  // NOLINT
    return def;
  }

  template <typename T = std::string>
  T Get(const Opt& opt) const {
    if (auto iter = values_.find(opt); iter != values_.end())
      return Conv<T>(iter->second);  // NOLINT
    error("missing option %s", opt.Desc().c_str());
    return {};
  }

  template <typename T = std::string>
  T Pos(std::size_t index) const {
    verify(index < pos_.size(), "index out of bounds");
    return Conv<T>(pos_[index]);
  }

  std::size_t PosSize() const { return pos_.size(); }

 private:
  // If the option is a flag, return the positive option first, and the inverted
  // one (--no-flag) second.
  std::pair<Opt, Opt> FlagPair(const Opt& opt) const;

  std::map<std::string, Opt> longname_;
  std::map<std::string, Opt> shortname_;
  std::set<Opt> opts_;
  std::map<Opt, std::string> values_;
  std::vector<std::string> pos_;
};

}  // namespace mrna

#endif  // UTIL_ARGPARSE_H_
