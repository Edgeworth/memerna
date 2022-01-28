// Copyright 2016 Eliot Courtney.
#include "util/argparse.h"

#include "util/string.h"

namespace mrna {

std::string Opt::Desc() const {
  std::string desc;
  if (!shortname_.empty() && !longname_.empty()) {
    desc += "-" + shortname_ + ", --" + longname_;
  } else if (!shortname_.empty()) {
    desc += "-" + shortname_;
  } else if (!longname_.empty()) {
    desc += "--" + longname_;
  }
  if (has_default_) desc += sfmt(" [%s]", default_arg_.c_str());
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

  verify(!opt.longname().empty() || !opt.shortname().empty(),
      "option must have either a shortname or a longname");

  if (!opt.longname().empty()) longname_.emplace(opt.longname(), opt);
  if (!opt.shortname().empty()) shortname_.emplace(opt.shortname(), opt);
  opts_.insert(opt);
  // Add default argument if necessary.
  if (auto iter = values_.find(opt); opt.has_default() && iter == values_.end())
    values_.emplace(opt, opt.default_arg());
}

std::string ArgParse::Usage() const {
  std::string usage = "Usage: \n";
  for (const auto& arg : opts_) usage += sfmt("  %s\n", arg.Desc().c_str());
  return usage;
}

std::string ArgParse::Parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    bool is_opt = arg[0] == '-';
    while (*arg == '-') ++arg;
    bool is_short = arg - argv[i] == 1;

    if (!is_opt) {
      positional_.push_back(arg);
    } else {
      auto& map = is_short ? shortname_ : longname_;
      auto iter = map.find(arg);
      if (iter == map.end()) return sfmt("unknown option %s", argv[i]);

      const auto& opt = iter->second;
      if (opt.has_arg()) {
        if (i + 1 == argc) return sfmt("missing argument for option %s", opt.Desc().c_str());
        values_[opt] = argv[++i];
      } else {
        values_[opt] = arg;
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
    fprintf(stderr, "%s\n%s\n", ret.c_str(), Usage().c_str());
    exit(1);
  }
}

const Opt& ArgParse::Lookup(const std::string& name) const {
  if (auto iter = longname_.find(name); iter != longname_.end()) return iter->second;  // NOLINT
  if (auto iter = shortname_.find(name); iter != shortname_.end()) return iter->second;  // NOLINT
  error("unregistered argument %s", name.c_str());
}

bool ArgParse::Has(const Opt& opt) const { return values_.contains(opt); }

std::string ArgParse::Get(const Opt& opt) const {
  if (auto iter = values_.find(opt); iter != values_.end()) return iter->second;  // NOLINT
  error("missing option %s", opt.Desc().c_str());
  return "";
}

}  // namespace mrna
