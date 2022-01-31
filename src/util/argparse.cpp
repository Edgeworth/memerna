// Copyright 2016 Eliot Courtney.
#include "util/argparse.h"

#include <cstdio>
#include <cstdlib>

#include "util/string.h"

namespace mrna {

std::string Opt::Desc() const {
  std::string desc;
  auto longname = longname_;

  // Handle flag --no.
  if (!longname.empty() && kind_ == FLAG) longname = longname + "/--no-" + longname;

  if (!shortname_.empty() && !longname.empty()) {
    desc += "-" + shortname_ + ", --" + longname;
  } else if (!shortname_.empty()) {
    desc += "-" + shortname_;
  } else if (!longname.empty()) {
    desc += "--" + longname;
  }
  if (has_default_) desc += sfmt(" [%s]", default_.c_str());
  if (!choices_.empty()) {
    desc += " (";
    for (auto iter = choices_.begin(); iter != choices_.end(); ++iter) {
      if (iter != choices_.begin()) desc += ", ";
      desc += *iter;
    }
    desc += ")";
  }
  if (!has_default_ && !choices_.empty()) desc += " <arg>";
  if (!help_.empty()) desc += ": " + help_;
  return desc;
}

void ArgParse::RegisterOpt(const Opt& opt) {
  // Check if this conflicts with existing arguments.
  if (auto iter = longname_.find(opt.longname()); iter != longname_.end())
    verify(opt == iter->second, "conflicting option registered with longname %s",
        opt.longname().c_str());
  if (auto iter = shortname_.find(opt.shortname()); iter != shortname_.end())
    verify(opt == iter->second, "conflicting option registered with shortname %s",
        opt.shortname().c_str());

  // If opt is a flag and has a longname, add the inversion to longname_ map.
  bool has_inversion = !opt.longname().empty() && opt.kind() == Opt::FLAG;
  std::string inverted_longname = "no-" + opt.longname();
  if (auto iter = longname_.find(inverted_longname); has_inversion && iter != longname_.end())
    verify(opt == iter->second, "conflicting option registered with longname %s",
        inverted_longname.c_str());

  verify(!opt.longname().empty() || !opt.shortname().empty(),
      "option must have either a shortname or a longname");

  if (!opt.longname().empty()) longname_.emplace(opt.longname(), opt);
  if (has_inversion)
    longname_.emplace(inverted_longname, Opt(opt).LongName(inverted_longname).Hidden());
  if (!opt.shortname().empty()) shortname_.emplace(opt.shortname(), opt);
  opts_.insert(opt);
  // Add default argument if necessary.
  if (auto iter = values_.find(opt); opt.has_default() && iter == values_.end())
    values_.emplace(opt, opt.default_arg());
}

std::string ArgParse::Usage() const {
  std::string usage = "Usage: \n";
  for (const auto& opt : opts_) {
    if (opt.hidden()) continue;
    usage += sfmt("  %s\n", opt.Desc().c_str());
  }
  return usage;
}

std::string ArgParse::Parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    const char* s = argv[i];
    bool is_opt = s[0] == '-';
    while (*s == '-') ++s;
    bool is_short = s - argv[i] == 1;

    if (!is_opt) {
      pos_.push_back(s);
    } else {
      auto& map = is_short ? shortname_ : longname_;
      auto iter = map.find(s);
      if (iter == map.end()) return sfmt("unknown option %s", argv[i]);

      const auto& opt = iter->second;
      if (opt.kind() == Opt::ARG) {
        if (i + 1 == argc) return sfmt("missing argument for option %s", opt.Desc().c_str());
        values_[opt] = argv[++i];
      } else {
        auto pair = FlagPair(opt);
        bool on = !opt.IsInverted();
        values_[pair.first] = convert(on);
        values_[pair.second] = convert(!on);
        printf("%s %s\n", values_[pair.first].c_str(), values_[pair.second].c_str());
      }
    }
  }
  for (const auto& opt : opts_) {
    bool has = Has(opt);
    if (opt.required() && !has) return sfmt("missing required option %s", opt.Desc().c_str());
    if (has && !opt.choices().empty() && !opt.choices().contains(Get(opt)))
      return sfmt("unrecognised argument for option %s", opt.Desc().c_str());
  }
  return "";
}

void ArgParse::ParseOrExit(int argc, char** argv) {
  const auto ret = Parse(argc, argv);
  if (!ret.empty()) {
    std::fprintf(stderr, "%s\n%s\n", ret.c_str(), Usage().c_str());
    std::exit(1);
  }
}

const Opt& ArgParse::Lookup(const std::string& name) const {
  if (auto iter = longname_.find(name); iter != longname_.end()) return iter->second;  // NOLINT
  if (auto iter = shortname_.find(name); iter != shortname_.end()) return iter->second;  // NOLINT
  error("unregistered argument %s", name.c_str());
}

bool ArgParse::Has(const Opt& opt) const { return values_.contains(opt); }

std::pair<Opt, Opt> ArgParse::FlagPair(const Opt& opt) const {
  verify(opt.kind() == Opt::FLAG, "option %s is not a flag", opt.longname().c_str());
  if (opt.IsInverted()) return {Lookup(opt.longname().substr(3)), opt};
  return {opt, Lookup("no-" + opt.longname())};
}

}  // namespace mrna
