// Copyright 2016 Eliot Courtney.
#ifndef UTIL_ARGPARSE_H_
#define UTIL_ARGPARSE_H_

#include <cstddef>
#include <map>
#include <optional>
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

  constexpr Opt& LongName(std::string n) {
    longname_ = std::move(n);
    return *this;
  }

  constexpr Opt& ShortName(std::string n) {
    shortname_ = std::move(n);
    return *this;
  }

  Opt& Choice(std::set<std::string> c) {
    verify(kind_ == ARG, "cannot set choices for flag");
    choices_ = std::move(c);
    return *this;
  }

  template <typename T>
  Opt& ChoiceEnum() {
    verify(kind_ == ARG, "cannot set choices for flag");
    auto choices = EnumNames<T>();
    choices_ = std::set(choices.begin(), choices.end());
    return *this;
  }

  template <typename T>
  Opt& Default(const T& d) {
    has_default_ = true;
    default_ = Conv(d);
    return *this;
  }

  Opt& AllChoicesAsDefault() {
    verify(!choices_.empty(), "choices must be set already");
    verify(multiple_, "multiple must be set");
    has_default_ = true;
    default_ = Join(choices_, ",");
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

  constexpr Opt& Multiple() {
    multiple_ = true;
    return *this;
  }

  Opt& Help(std::string h) {
    help_ = std::move(h);
    return *this;
  }

  [[nodiscard]] std::string Desc() const;

  [[nodiscard]] bool IsInverted() const { return longname_.rfind("no-", 0) == 0; }

  auto operator<=>(const Opt&) const = default;

  [[nodiscard]] Kind kind() const { return kind_; }
  [[nodiscard]] const std::string& longname() const { return longname_; }
  [[nodiscard]] const std::string& shortname() const { return shortname_; }
  [[nodiscard]] const std::string& help() const { return help_; }
  [[nodiscard]] const std::string& default_arg() const { return default_; }
  [[nodiscard]] const std::set<std::string>& choices() const { return choices_; }

  [[nodiscard]] bool has_default() const { return has_default_; }
  [[nodiscard]] bool required() const { return required_; }
  [[nodiscard]] bool hidden() const { return hidden_; }
  [[nodiscard]] bool multiple() const { return multiple_; }

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
  bool multiple_ = false;
};

class ArgParse {
 public:
  ArgParse() = default;
  ~ArgParse() = default;
  ArgParse(const ArgParse&) = delete;
  ArgParse& operator=(const ArgParse&) = delete;

  ArgParse(ArgParse&&) = default;
  ArgParse& operator=(ArgParse&&) = default;

  void RegisterOpt(const Opt& opt);
  [[nodiscard]] std::string Usage() const;

  std::string Parse(int argc, char* argv[]);
  void ParseOrExit(int argc, char* argv[]);

  [[nodiscard]] const Opt& Lookup(const std::string& name) const;
  [[nodiscard]] const std::vector<std::string>& Pos() const { return pos_; }

  // Note that flags set to false via --no-flag will have Has return true.
  // It returns whether or not this option was specified. Use GetOr instead
  // for flags.
  [[nodiscard]] bool Has(const Opt& opt) const;

  // Prefer this for setting variables over Get<T>(), if those variables
  // have default values.
  template <typename T>
  void MaybeSet(const Opt& opt, T* val) const {
    if (auto iter = values_.find(opt); iter != values_.end())
      *val = Conv<T>(iter->second);  // NOLINT
  }

  template <typename T>
  [[nodiscard]] std::optional<T> MaybeGet(const Opt& opt) const {
    if (auto iter = values_.find(opt); iter != values_.end())
      return Conv<T>(iter->second);  // NOLINT
    return std::nullopt;
  }

  // Useful mainly with flags where not specifying them means false (or whatever
  // you pass as the default).
  template <typename T = bool>
  [[nodiscard]] T GetOr(const Opt& opt, T def = T()) const {
    if (auto iter = values_.find(opt); iter != values_.end())
      return Conv<T>(iter->second);  // NOLINT
    return def;
  }

  template <typename T = std::string>
  [[nodiscard]] T Get(const Opt& opt) const {
    if (auto iter = values_.find(opt); iter != values_.end())
      return Conv<T>(iter->second);  // NOLINT
    fatal("missing option {}", opt.Desc());
    return {};
  }

  template <typename T = std::string>
  [[nodiscard]] std::vector<T> GetMultiple(const Opt& opt) const {
    auto iter = values_.find(opt);
    if (iter == values_.end()) fatal("missing option {}", opt.Desc());
    std::vector<T> values;
    for (const auto& s : Split(iter->second, ",")) values.emplace_back(Conv<T>(s));
    return values;
  }

  template <typename T = std::string>
  [[nodiscard]] std::vector<T> GetMultipleOr(
      const Opt& opt, std::vector<T> def = std::vector<T>()) const {
    if (!values_.contains(opt)) return def;
    return GetMultiple<T>(opt);
  }

  template <typename T = std::string>
  [[nodiscard]] T Pos(std::size_t index) const {
    verify(index < pos_.size(), "index out of bounds");
    return Conv<T>(pos_[index]);
  }

  [[nodiscard]] std::size_t PosSize() const { return pos_.size(); }

 private:
  // If the option is a flag, return the positive option first, and the inverted
  // one (--no-flag) second.
  [[nodiscard]] std::pair<Opt, Opt> FlagPair(const Opt& opt) const;

  std::map<std::string, Opt> longname_;
  std::map<std::string, Opt> shortname_;
  std::vector<Opt> opts_;  // In order of registration.
  std::map<Opt, std::string> values_;
  std::vector<std::string> pos_;
};

}  // namespace mrna

#endif  // UTIL_ARGPARSE_H_
