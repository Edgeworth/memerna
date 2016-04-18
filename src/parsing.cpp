#include "parsing.h"

namespace memerna {
namespace parsing {

rna_t ParseRnaFromString(const std::string& s) {
  rna_t rna(s.size());
  for (int i = 0; i < int(s.size()); ++i) {
    rna[i] = CharToBase(s[i]);
  }
  return rna;
}

folded_rna_t ParseViennaRna(const std::string& rna_str, const std::string& pairs_str) {
  assert(rna_str.size() == pairs_str.size());
  rna_t rna = ParseRnaFromString(rna_str);
  std::vector<int> pairs(rna_str.size(), -1);
  std::stack<int> s;
  for (int i = 0; i < int(rna_str.size()); ++i) {
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


void Parse2x2FromFile(const std::string& filename, energy::energy_t (& output)[4][4][4][4]) {
  FILE* fp = fopen(filename.c_str(), "r");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    assert(a != -1 && b != -1 && c != -1 && d != -1);
    int res = fscanf(fp, " %d ", &output[a][b][c][d]);
    assert(res == 1);
  }
  fclose(fp);
}

void ParseMapFromFile(const std::string& filename, std::unordered_map<std::string, energy::energy_t>& output) {
  FILE* fp = fopen(filename.c_str(), "r");
  char buf[1024];
  energy::energy_t energy;
  while (fscanf(fp, " %s %d ", buf, &energy) == 2)
    output[buf] = energy;
  fclose(fp);
}

void ParseInitiationEnergyFromFile(const std::string& filename, energy::energy_t (& output)[INITIATION_CACHE_SZ]) {
  FILE* fp = fopen(filename.c_str(), "r");
  energy::energy_t energy;
  int idx;
  while (fscanf(fp, "%d %d ", &idx, &energy) == 2)
    output[idx] = energy;
  fclose(fp);
}

void ParseHairpinMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(
      fp, "%d %d %d %d %d %d",
      &hairpin_uu_ga_first_mismatch, &hairpin_gg_first_mismatch,
      &hairpin_special_gu_closure, &hairpin_c3_loop, &hairpin_all_c_a, &hairpin_all_c_b);
  assert(res == 6);
  fclose(fp);
}

void ParseBulgeMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(fp, "%d", &bulge_special_c);
  assert(res == 1);
  fclose(fp);
}

void ParseInternalLoop1x1FromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    base_t e = CharToBase((char) fgetc(fp));
    base_t f = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    assert(a != -1 && b != -1 && c != -1 && d != -1 && e != -1 && f != -1);
    int res = fscanf(fp, " %d ", &internal_1x1[a][b][c][d][e][f]);
    assert(res == 1);
  }
  fclose(fp);
}

void ParseInternalLoop1x2FromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    base_t e = CharToBase((char) fgetc(fp));
    base_t f = CharToBase((char) fgetc(fp));
    base_t g = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    assert(a != -1 && b != -1 && c != -1 && d != -1 && e != -1 && f != -1 && g != -1);
    int res = fscanf(fp, " %d ", &internal_1x2[a][b][c][d][e][f][g]);
    assert(res == 1);
  }
  fclose(fp);
}

void ParseInternalLoop2x2FromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    base_t e = CharToBase((char) fgetc(fp));
    base_t f = CharToBase((char) fgetc(fp));
    base_t g = CharToBase((char) fgetc(fp));
    base_t h = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    assert(a != -1 && b != -1 && c != -1 && d != -1 && e != -1 && f != -1 && g != -1 && h != -1);
    int res = fscanf(fp, " %d ", &internal_2x2[a][b][c][d][e][f][g][h]);
    assert(res == 1);
  }
  fclose(fp);
}

void ParseInternalLoopMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(fp, "%d %d %d", &internal_asym, &internal_augu_penalty, &internal_mismatch_1xk);
  assert(res == 3);
  fclose(fp);
}

void ParseMultiloopT99MiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(fp, "%d %d %d", &multiloop_t99_a, &multiloop_t99_b, &multiloop_t99_c);
  assert(res == 3);
  fclose(fp);
}

void ParseMultiloopHackMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(fp, "%d %d", &multiloop_hack_a, &multiloop_hack_b);
  assert(res == 2);
  fclose(fp);
}

void ParseDangleDataFromFile(const std::string& filename, energy::energy_t (& output)[4][4][4]) {
  FILE* fp = fopen(filename.c_str(), "r");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    assert(a != -1 && b != -1 && c != -1);
    int res = fscanf(fp, " %d ", &output[a][b][c]);
    assert(res == 1);
  }
  fclose(fp);
}

void ParseCoaxialMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(fp, "%d %d %d", &coax_mismatch_non_contiguous, &coax_mismatch_wc_bonus, &coax_mismatch_gu_bonus);
  assert(res == 3);
  fclose(fp);
}

}
}
