#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include <stack>
#include <cassert>

namespace memerna {
namespace base {

typedef int8_t base_t;
typedef int32_t energy_t;
typedef std::vector<base_t> rna_t;

struct folded_rna_t {
    rna_t rna;
    std::vector<int> pairs;
};

namespace internal {

energy_t ComputeEnergyInternal(const folded_rna_t& folded_rna, int st, int en) {
  return 0;
}

}

base_t BaseFromChar(char c) {
  switch (c) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'U':
      return 3;
    default:
      return -1;
  }
}

char BaseToChar(base_t b) {
  if (b < 0 || b > 4) return '?';
  return "ACGU"[b];
}

energy_t ComputeEnergy(const folded_rna_t& folded_rna) {
  energy_t energy = 0;
  for (auto i :  folded_rna.pairs) {
    energy += i != -1;
  }
  assert(energy % 2 == 0);
  return energy / 2;
}

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
}
}

using namespace memerna::base;

int main(int argc, char* argv[]) {
  printf("meme");
  assert(argc == 3);
  folded_rna_t frna = ParseViennaRna(argv[1], argv[2]);
  printf("Computed energy: %d", ComputeEnergy(frna));
}
