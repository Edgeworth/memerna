#include "parsing.h"
#include "globals.h"

namespace memerna {
namespace parsing {

rna_t ParseRnaFromString(const std::string& s) {
  rna_t rna(s.size());
  for (int i = 0; i < s.size(); ++i) {
    rna[i] = BaseFromChar(s[i]);
  }
  return rna;
}

folded_rna_t ParseViennaRna(const std::string& rna_str, const std::string& pairs_str) {
  rna_t rna = ParseRnaFromString(rna_str);
  std::vector<int> pairs(rna_str.size(), -1);
  std::stack<int> s;
  for (int i = 0; i < rna_str.size(); ++i) {
    if (pairs_str[i] == '(') {
      s.push(i);
    } else if (pairs_str[i] == ')') {
      assert(!s.empty());
      pairs[i] = s.top();
      pairs[s.top()] = i;
      s.pop();
    }
  }
  return {rna, pairs};
}


void ParseStackingEnergiesFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  for (int i = 0; i < 256; ++i) {
    base_t a = BaseFromChar((char) fgetc(fp));
    base_t b = BaseFromChar((char) fgetc(fp));
    base_t c = BaseFromChar((char) fgetc(fp));
    base_t d = BaseFromChar((char) fgetc(fp));
    assert(a != -1 && b != -1 && c != -1 && d != -1);
    int res = fscanf(fp, " %d ", &stacking_e[a][b][c][d]);
    assert(res == 1);
  }
  fclose(fp);
}

}
}
