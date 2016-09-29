#include <set>
#include "common.h"
#include "argparse.h"
#include "splaymap.h"

using namespace memerna;

int main(int argc, char** argv) {
  ArgParse argparse({
      {"h", {"help"}},
      {"afl", {"afl mode"}},
      {"r", ArgParse::option_t("load range from 1 until r").Arg()}
  });
  argparse.ParseOrExit(argc, argv);
  if (argparse.HasFlag("h")) {
    printf("Commands:\n Insert: i <val>\n Delete: d <val>\n "
        "Search: s <val>\n Range: a <min> <max>\n");
    printf("%s\n", argparse.Usage().c_str());
    return 0;
  }
  const bool afl = argparse.HasFlag("afl");
#ifdef __AFL_HAVE_MANUAL_CONTROL
  __AFL_INIT();
    while (__AFL_LOOP(100000)) {
#endif

  SplayMap<int, int> h;
  std::set<int> s;
  if (argparse.HasFlag("r")) {
    int r = atoi(argparse.GetOption("r").c_str());
    for (int i = 1; i <= r; ++i) {
      h.Insert(i, i);
      s.insert(i);
    }
    printf("%s\n", h.Describe().c_str());
  }

  while(1) {
    if (!afl) printf("> ");
    char op = 0;
    int val = 0;
    int res = scanf(" %c %d", &op, &val);
    if (res >= 0) {
      if (res != 2) {
        if (!afl) printf("Invalid input.\n");
        continue;
      }
    } else {
      break;
    }

    bool success = false;
    bool set_success = false;
    switch (op) {
      case 'i':
        success = h.Insert(val, val);
        set_success = (s.count(val) == 0);
        s.insert(val);
        break;
      case 'd':
        success = h.Delete(val);
        set_success = (s.count(val) > 0);
        s.erase(val);
        break;
      case 's':
        success = h.Find(val);
        set_success = (s.find(val) != s.end());
        break;
      default:
        if (!afl) printf("Invalid input.\n");
        continue;
    }
    if (!afl) {
      if (success) printf("Operation successful.\n");
      else printf("Operation failed.\n");
      printf("%s\n", h.Describe().c_str());
    }
    verify_expr(success == set_success, "bug");
    verify_expr(s.size() == h.Size(), "bug");
    auto keys = h.Keys();
    verify_expr(std::equal(s.begin(), s.end(), keys.begin(), keys.end()), "bug");
  }
#ifdef __AFL_HAVE_MANUAL_CONTROL
  }
#endif
}
