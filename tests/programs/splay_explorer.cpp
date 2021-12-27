// Copyright 2016 E.
#include <iomanip>
#include <set>
#include <sstream>

#include "util/argparse.h"
#include "util/splaymap.h"

enum class OpResult { INVALID, FAILURE, SUCCESS };

OpResult DoOperation(char op, int val, mrna::SplayMap<int, int>& h, std::set<int>& s) {
  bool splay_success = false, set_success = false, invalid = false;
  switch (op) {
  case 'i':
    splay_success = h.Insert(val, val);
    set_success = (s.count(val) == 0);
    s.insert(val);
    break;
  case 'd':
    splay_success = h.Delete(val);
    set_success = (s.count(val) > 0);
    s.erase(val);
    break;
  case 's':
    splay_success = h.Find(val);
    set_success = (s.find(val) != s.end());
    break;
  default: invalid = true;
  }
  auto keys = h.Keys();
  if (splay_success != set_success || s.size() != h.Size() ||
      !std::equal(s.begin(), s.end(), keys.begin(), keys.end()))
    std::abort();
  if (invalid) return OpResult::INVALID;
  if (splay_success) return OpResult::SUCCESS;
  return OpResult::FAILURE;
}

void DoAfl() {
#ifdef __AFL_HAVE_MANUAL_CONTROL
  __AFL_INIT();
  while (__AFL_LOOP(100000)) {
#endif
    mrna::SplayMap<int, int> h;
    std::set<int> s;

    std::string data;
    std::size_t len;
    char buf[4096];
    while ((len = fread(buf, 1, sizeof(buf), stdin)) > 0) data += std::string(buf, len);
    std::stringstream ss(data);
    char op = 0;
    int val = 0;
    while (ss >> std::ws >> op >> std::ws >> val) DoOperation(op, val, h, s);
    printf("%s\n", h.Describe().c_str());
#ifdef __AFL_HAVE_MANUAL_CONTROL
  }
#endif
}

void DoInteractive(int r) {
  mrna::SplayMap<int, int> h;
  std::set<int> s;

  for (int i = 1; i <= r; ++i) {
    h.Insert(i, i);
    s.insert(i);
  }

  char op = 0;
  int val = 0;
  while (1) {
    printf("%s\n> ", h.Describe().c_str());
    int res = scanf(" %c %d", &op, &val);
    if (res < 0) break;
    if (res == 2) {
      auto op_res = DoOperation(op, val, h, s);
      if (op_res == OpResult::INVALID)
        printf("Invalid input\n");
      else if (op_res == OpResult::SUCCESS)
        printf("Operation successful.\n");
      else
        printf("Operation failed.\n");
    } else {
      printf("Invalid input.\n");
    }
  }
}

int main(int argc, char** argv) {
  mrna::ArgParse args({{"h", {"help"}}, {"afl", {"afl mode"}},
      {"r", mrna::opt_t("load range from 1 until r").Arg("-1")}});
  args.ParseOrExit(argc, argv);
  if (args.HasFlag("h")) {
    printf("Commands:\n Insert: i <val>\n Delete: d <val>\n "
           "Search: s <val>\n Range: a <min> <max>\n");
    printf("%s\n", args.Usage().c_str());
    return 0;
  }
  const bool afl = args.HasFlag("afl");
  const int r = atoi(args.GetOption("r").c_str());
  verify(!(afl && (r != -1)), "incompatible options");

  if (afl)
    DoAfl();
  else
    DoInteractive(r);
}
