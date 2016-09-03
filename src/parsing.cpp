#include "parsing.h"

namespace memerna {
namespace parsing {

using namespace energy;

primary_t StringToPrimary(const std::string& s) {
  primary_t r(s.size());
  for (int i = 0; i < int(s.size()); ++i) {
    r[i] = CharToBase(s[i]);
    verify_expr(r[i] != -1, "unexpected base %c", s[i]);
  }
  return r;
}

secondary_t ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str) {
  verify_expr(prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
  return {StringToPrimary(prim_str), DotBracketToPairs(pairs_str)};
}

std::vector<int> DotBracketToPairs(const std::string& pairs_str) {
  std::vector<int> pairs(pairs_str.size(), -1);
  std::stack<int> s;
  for (int i = 0; i < int(pairs_str.size()); ++i) {
    if (pairs_str[i] == '(') {
      s.push(i);
    } else if (pairs_str[i] == ')') {
      verify_expr(!s.empty(), "unmatched bracket");
      pairs[i] = s.top();
      pairs[s.top()] = i;
      s.pop();
    }
  }
  return pairs;
}

std::string PairsToDotBracket(const std::vector<int>& pairs) {
  std::string s(pairs.size(), '.');
  for (int i = 0; i < int(pairs.size()); ++i) {
    if (pairs[i] == -1) continue;
    if (pairs[i] < i) s[i] = ')';
    else s[i] = '(';
  }
  return s;
}

std::string ComputedToCtdString(const computed_t& computed) {
  const auto& p = computed.s.p;
  std::string s = PairsToDotBracket(p);
  for (int i = 0; i < int(p.size()); ++i) {
    switch (computed.base_ctds[i]) {
      case CTD_NA:
      case CTD_UNUSED:
        break;
      case CTD_3_DANGLE:
        s[p[i] + 1] = '3';
        break;
      case CTD_5_DANGLE:
        s[i - 1] = '5';
        break;
      case CTD_TERMINAL_MISMATCH:
        s[i - 1] = 'm';
        s[p[i] + 1] = 'm';
        break;
      case CTD_LEFT_MISMATCH_COAX_WITH_NEXT:
        s[i] = '>';
        s[i - 1] = 'm';
        s[p[i] - 1] = 'm';
        break;
      case CTD_LEFT_MISMATCH_COAX_WITH_PREV:
        s[i] = '<';
        break;
      case CTD_RIGHT_MISMATCH_COAX_WITH_NEXT:
        s[i] = '>';
        break;
      case CTD_RIGHT_MISMATCH_COAX_WITH_PREV:
        s[i] = '<';
        s[i - 1] = 'm';
        s[p[i] - 1] = 'm';
        break;
      case CTD_FLUSH_COAX_WITH_NEXT:
        s[i] = '>';
        break;
      case CTD_FLUSH_COAX_WITH_PREV:
        s[i] = '<';
        break;
      default:
        verify_expr(false, "bug");
    }
  }
  return s;
}

std::string PrimaryToString(const primary_t& r) {
  std::string s;
  s.resize(r.size());
  for (int i = 0; i < int(r.size()); ++i) {
    s[i] = BaseToChar(r[i]);
  }
  return s;
}

void Parse2x2FromFile(const std::string& filename, energy_t (& output)[4][4][4][4]) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    verify_expr(a != -1 && b != -1 && c != -1 && d != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &output[a][b][c][d]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseMapFromFile(const std::string& filename, std::unordered_map<std::string, energy_t>& output) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  char buf[1024];
  energy_t energy;
  while (fscanf(fp, " %s %d ", buf, &energy) == 2)
    output[buf] = energy;
  fclose(fp);
}

void ParseInitiationEnergyFromFile(const std::string& filename, energy_t (& output)[INITIATION_CACHE_SZ]) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  energy_t energy;
  int idx;
  while (fscanf(fp, "%d %d ", &idx, &energy) == 2)
    output[idx] = energy;
  fclose(fp);
}

void ParseInternalLoop1x1FromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    base_t e = CharToBase((char) fgetc(fp));
    base_t f = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    verify_expr(a != -1 && b != -1 && c != -1 && d != -1 && e != -1 && f != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &g_internal_1x1[a][b][c][d][e][f]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseInternalLoop1x2FromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    base_t e = CharToBase((char) fgetc(fp));
    base_t f = CharToBase((char) fgetc(fp));
    base_t g = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    verify_expr(a != -1 && b != -1 && c != -1 && d != -1 &&
        e != -1 && f != -1 && g != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &g_internal_1x2[a][b][c][d][e][f][g]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseInternalLoop2x2FromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
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
    verify_expr(a != -1 && b != -1 && c != -1 && d != -1 && e != -1 &&
        f != -1 && g != -1 && h != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &g_internal_2x2[a][b][c][d][e][f][g][h]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseDangleDataFromFile(const std::string& filename, energy_t (& output)[4][4][4]) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    if (a == -1) break;
    verify_expr(a != -1 && b != -1 && c != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &output[a][b][c]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");

#define READ_DATA(var) \
  do { \
    while (1) { \
      std::string line = sgetline(fp); \
      verify_expr(line.size() > 0, "unexpected EOF or error"); \
      if (line[0] == '/' || line[0] == '\n') continue; \
      verify_expr(sscanf(line.c_str(), "%d", &var) == 1, "expected int"); \
      break; \
    } \
  } while (0)


  // Bulge loops.
  READ_DATA(g_bulge_special_c);

  // Coaxial stacking.
  READ_DATA(g_coax_mismatch_non_contiguous);
  READ_DATA(g_coax_mismatch_wc_bonus);
  READ_DATA(g_coax_mismatch_gu_bonus);

  // Hairpin loops.
  READ_DATA(g_hairpin_uu_ga_first_mismatch);
  READ_DATA(g_hairpin_gg_first_mismatch);
  READ_DATA(g_hairpin_special_gu_closure);
  READ_DATA(g_hairpin_c3_loop);
  READ_DATA(g_hairpin_all_c_a);
  READ_DATA(g_hairpin_all_c_b);

  // Internal loops.
  READ_DATA(g_internal_asym);
  READ_DATA(g_internal_augu_penalty);

  // Multiloop data.
  READ_DATA(g_multiloop_hack_a);
  READ_DATA(g_multiloop_hack_b);

  // AU/GU penalty
  READ_DATA(g_augu_penalty);
#undef READ_DATA

  fclose(fp);
}

}
}
