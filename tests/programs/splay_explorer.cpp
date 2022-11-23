// Copyright 2016 Eliot Courtney.
#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "options.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/splaymap.h"

#ifdef __AFL_FUZZ_TESTCASE_LEN
#include <unistd.h>  // For __AFL_FUZZ_TESTCASE_LEN

__AFL_FUZZ_INIT();
#endif

enum class OpResult { INVALID, FAILURE, SUCCESS };

OpResult DoOperation(char op, int val, mrna::SplayMap<int, int>* h, std::set<int>* s) {
  bool splay_success = false;
  bool set_success = false;
  bool invalid = false;
  switch (op) {
  case 'i':
    splay_success = h->Insert(val, val);
    set_success = (s->count(val) == 0);
    s->insert(val);
    break;
  case 'd':
    splay_success = h->Delete(val);
    set_success = (s->count(val) > 0);
    s->erase(val);
    break;
  case 's':
    splay_success = h->Find(val);
    set_success = (s->find(val) != s->end());
    break;
  default: invalid = true;
  }
  auto keys = h->Keys();
  if (splay_success != set_success || s->size() != h->Size() ||
      !std::equal(s->begin(), s->end(), keys.begin(), keys.end()))
    std::abort();
  if (invalid) return OpResult::INVALID;
  if (splay_success) return OpResult::SUCCESS;
  return OpResult::FAILURE;
}

void DoAfl() {
#ifdef __AFL_FUZZ_TESTCASE_LEN
  __AFL_INIT();
  // This must be after __AFL_INIT and before __AFL_LOOP.
  auto buf = reinterpret_cast<const char*>(__AFL_FUZZ_TESTCASE_BUF);
  while (__AFL_LOOP(100000)) {
    int len = __AFL_FUZZ_TESTCASE_LEN;
    std::string data(buf, len);

    mrna::SplayMap<int, int> h;
    std::set<int> s;
    std::stringstream ss(data);
    char op = 0;
    int val = 0;
    while (ss >> std::ws >> op >> std::ws >> val) DoOperation(op, val, &h, &s);
    std::cout << h.Describe() << '\n';
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
  while (true) {
    std::cout << h.Describe() << "\n>";
    if (!(std::cin >> op >> val)) {
      std::cout << "Invalid input.\n";
      break;
    }
    auto op_res = DoOperation(op, val, &h, &s);
    if (op_res == OpResult::INVALID)
      std::cout << "Invalid input\n";
    else if (op_res == OpResult::SUCCESS)
      std::cout << "Operation successful.\n";
    else
      std::cout << "Operation failed.\n";
  }
}

inline const auto OPT_RANGE =
    mrna::Opt(mrna::Opt::ARG).ShortName("r").Default("-1").Help("load range from 1 until r");

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  args.RegisterOpt(mrna::OPT_AFL);
  args.RegisterOpt(OPT_RANGE);
  args.ParseOrExit(argc, argv);

  std::cout
      << "Commands:\n Insert: i <val>\n Delete: d <val>\nSearch: s <val>\n Range: a <min> <max>\n";
  const bool afl = args.GetOr(mrna::OPT_AFL);
  const int r = args.Get<int>(OPT_RANGE);
  verify(!(afl && (r != -1)), "incompatible options");

  if (afl)
    DoAfl();
  else
    DoInteractive(r);
}
