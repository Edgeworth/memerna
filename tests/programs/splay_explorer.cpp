#include "common.h"
#include "argparse.h"
#include "splaymap.h"

using namespace memerna;

int main(int argc, char** argv) {
  ArgParse argparse({
      {"h", {"help"}},
      {"r", ArgParse::option_t("load range from 1 until r").Arg()}
  });
  argparse.ParseOrExit(argc, argv);
  if (argparse.HasFlag("h")) {
    printf("Commands:\n Insert: i <val>\n Delete: d <val>\n Search: s <val>\n");
    printf("%s\n", argparse.Usage().c_str());
    return 0;
  }

  SplayMap<int, int> h;
  if (argparse.HasFlag("r")) {
    int r = atoi(argparse.GetOption("r").c_str());
    for (int i = 1; i <= r; ++i)
      h.Insert(i, i);
    printf("%s\n", h.Describe().c_str());
  }

  while(1) {
    printf("> ");
    char op = 0;
    int val = 0;
    int res = scanf(" %c %d", &op, &val);
    if (res == 0) {
      printf("Invalid input.\n");
      continue;
    } else if (res != 2){
      break;
    }

    switch (op) {
      case 'i':
        if (h.Insert(val, val)) printf("Inserted.\n");
        else printf("Not inserted.\n");
        break;
      case 'd':
        if (h.Delete(val)) printf("Deleted.\n");
        else printf("Not deleted.\n");
        break;
      case 's':
        if (h.Find(val)) printf("Found.\n");
        else printf("Not found.\n");
        break;
      default:
        printf("Invalid input.\n");
        continue;
    }
    printf("%s\n", h.Describe().c_str());
  }
}
