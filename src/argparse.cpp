#include "argparse.h"

namespace memerna {

std::string ArgParse::option_t::Desc() const {
  auto res = desc;
  if (has_default)
    res += sfmt(" [%s]", default_arg.c_str());
  if (choices.size()) {
    res += " (";
    for (auto iter = choices.begin(); iter != choices.end(); ++iter) {
      if (iter != choices.begin())
        res += ", ";
      res += *iter;
    }
    res += ")";
  }
  return res;
}

std::string ArgParse::Parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    bool is_flag = arg[0] == '-';
    while (*arg == '-') ++arg;

    if (is_flag) {
      if (possible_args.count(arg) == 0)
        return sfmt("unknown argument %s", argv[i]);

      if (possible_args[arg].has_arg) {
        if (i + 1 == argc)
          return sfmt("missing argument for flag %s", arg);
        flags[arg] = argv[++i];
      } else {
        flags[arg] = arg;
      }
    } else {
      positional.push_back(std::move(arg));
    }
  }
  for (const auto& argpair : possible_args) {
    const auto& flag = argpair.first;
    const auto& arg = argpair.second;
    verify_expr(!arg.has_arg || arg.has_default, "bad option somehow");
    if (arg.has_arg && flags.count(flag) == 0 && !arg.has_default && arg.required)
      return sfmt("missing argument for flag %s", flag.c_str());
    if (!arg.choices.empty() && arg.choices.count(GetOption(flag)) == 0)
      return sfmt("unrecognised argument for flag %s", flag.c_str());
  }
  return "";
}

std::string ArgParse::Usage() const {
  std::string usage = "Flags: \n";
  for (const auto& arg : possible_args) {
    usage += sfmt("  -%s: %s\n", arg.first.c_str(), arg.second.Desc().c_str());
  }
  return usage;
}

void ArgParse::AddOptions(const std::map<std::string, ArgParse::option_t>& possible_args_) {
  for (const auto& argpair : possible_args_) {
    verify_expr(possible_args.count(argpair.first) == 0, "duplicate argument %s", argpair.first.c_str());
    possible_args.insert(argpair);
  }
}

void ArgParse::ParseOrExit(int argc, char** argv) {
  auto ret = Parse(argc, argv);
  verify_expr(ret.size() == 0, "%s\n%s\n", ret.c_str(), Usage().c_str());
}

fold::fold_fn_t* FoldFunctionFromArgParse(const ArgParse& argparse) {
  fold::fold_fn_t* fold_fn = nullptr;
  auto opt = argparse.GetOption("alg");
  if (opt == "0")
    fold_fn = &fold::Fold0;
  else if (opt == "1")
    fold_fn = &fold::Fold1;
  else if (opt == "2")
    fold_fn = &fold::Fold2;
  else if (opt == "3")
    fold_fn = &fold::Fold3;
  else if (opt == "brute")
    fold_fn = &fold::FoldBruteForce;
  else
    verify_expr(false, "unknown fold option");
  return fold_fn;
}


}
