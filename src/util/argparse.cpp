// Copyright 2016 Eliot Courtney.
#include "util/argparse.h"

#include <fmt/core.h>
#include <spdlog/common.h>

#include <algorithm>
#include <cstdlib>
#include <exception>

#include "api/options.h"
#include "spdlog/spdlog.h"
#include "util/string.h"

namespace mrna {

std::string Opt::Desc() const {
  std::string desc;
  auto ln = longname_;

  // Handle flag --no.
  if (!ln.empty() && kind_ == FLAG) ln = ln + "/--no-" + ln;

  if (!shortname_.empty() && !ln.empty()) {
    desc += "-" + shortname_ + ", --" + ln;
  } else if (!shortname_.empty()) {
    desc += "-" + shortname_;
  } else if (!ln.empty()) {
    desc += "--" + ln;
  }
  if (has_default_) desc += fmt::format(" [{}]", default_);
  if (!choices_.empty()) {
    desc += " (";
    for (auto iter = choices_.begin(); iter != choices_.end(); ++iter) {
      if (iter != choices_.begin()) desc += ", ";
      desc += *iter;
    }
    desc += ")";
  }
  if (multiple_) desc += " (multiple)";
  if (!has_default_ && !choices_.empty()) desc += " <arg>";
  if (!help_.empty()) desc += ": " + help_;
  return desc;
}

void ArgParse::RegisterOpt(const Opt& opt) {
  // Check if this conflicts with existing arguments.
  if (auto iter = longname_.find(opt.longname()); iter != longname_.end())
    verify(opt == iter->second, "conflicting option registered with longname {}", opt.longname());
  if (auto iter = shortname_.find(opt.shortname()); iter != shortname_.end())
    verify(opt == iter->second, "conflicting option registered with shortname {}", opt.shortname());

  // Can't have multiple of a flag:
  verify(opt.kind() != Opt::FLAG || !opt.multiple(), "flag {} cannot be multiple", opt.Desc());

  // If opt is a flag and has a longname, add the inversion to longname_ map.
  const bool has_inversion = !opt.longname().empty() && opt.kind() == Opt::FLAG;
  std::string inverted_longname = "no-" + opt.longname();
  auto inverted_opt = Opt(opt).LongName(inverted_longname).Hidden();
  if (auto iter = longname_.find(inverted_longname); has_inversion && iter != longname_.end())
    verify(inverted_opt == iter->second, "conflicting option registered with longname {}",
        inverted_longname);

  verify(!opt.longname().empty() || !opt.shortname().empty(),
      "option must have either a shortname or a longname");

  if (!opt.longname().empty()) longname_.emplace(opt.longname(), opt);
  if (has_inversion) longname_.emplace(inverted_longname, inverted_opt);
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
    usage += fmt::format("  {}\n", opt.Desc());
  }
  return usage;
}

std::string ArgParse::Parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    const char* s = argv[i];
    const bool is_opt = s[0] == '-';
    while (*s == '-') ++s;
    const bool is_short = s - argv[i] == 1;

    if (!is_opt) {
      pos_.emplace_back(s);
    } else {
      auto& map = is_short ? shortname_ : longname_;
      auto iter = map.find(s);
      if (iter == map.end()) return fmt::format("unknown option {}", argv[i]);

      const auto& opt = iter->second;
      if (opt.kind() == Opt::ARG) {
        if (i + 1 == argc) return fmt::format("missing argument for option {}", opt.Desc());
        if (opt.multiple()) {
          values_[opt] = std::string(argv[++i]);
        } else {
          values_[opt] = argv[++i];
        }
      } else {
        auto pair = FlagPair(opt);
        const bool on = !opt.IsInverted();
        values_[pair.first] = Conv(on);
        values_[pair.second] = Conv(!on);
      }
    }
  }
  for (const auto& opt : opts_) {
    const bool has = Has(opt);
    if (opt.required() && !has) return fmt::format("missing required option {}", opt.Desc());
    if (has && !opt.choices().empty()) {
      std::vector<std::string> string_values;
      if (opt.multiple()) {
        string_values = GetMultiple(opt);
      } else {
        string_values.emplace_back(Get(opt));
      }
      for (const auto& v : string_values) {
        if (!opt.choices().contains(v))
          return fmt::format("unrecognised argument for option {}", opt.Desc());
      }
    }
  }
  return "";
}

void ArgParse::ParseOrExit(int argc, char** argv) {
  RegisterOpt(OPT_VERBOSE);
  const auto ret = Parse(argc, argv);
  if (!ret.empty()) {
    spdlog::critical("{}\n{}\n", ret, Usage());
    std::exit(1);  // NOLINT
  }

  if (Has(OPT_VERBOSE)) {
    spdlog::set_level(spdlog::level::debug);
  }
}

const Opt& ArgParse::Lookup(const std::string& name) const {
  if (auto iter = longname_.find(name); iter != longname_.end()) return iter->second;  // NOLINT
  if (auto iter = shortname_.find(name); iter != shortname_.end()) return iter->second;  // NOLINT
  fatal("unregistered argument {}", name);
}

bool ArgParse::Has(const Opt& opt) const { return values_.contains(opt); }

std::pair<Opt, Opt> ArgParse::FlagPair(const Opt& opt) const {
  verify(opt.kind() == Opt::FLAG, "option {} is not a flag", opt.longname());
  if (opt.IsInverted()) return {Lookup(opt.longname().substr(3)), opt};
  return {opt, Lookup("no-" + opt.longname())};
}

}  // namespace mrna
