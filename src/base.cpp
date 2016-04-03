#include "base.h"
#include "parsing.h"

namespace memerna {

void Init() {
  parsing::Parse2x2FromFile("data/stacking.data", stacking_e);
  parsing::Parse2x2FromFile("data/terminal.data", terminal_e);
  parsing::ParseMapFromFile("data/hairpin.data", hairpin_e);
  parsing::ParseInitiationEnergyFromFile("data/internal_initiation.data", internal_init);
  parsing::ParseInitiationEnergyFromFile("data/bulge_initiation.data", bulge_init);
  parsing::ParseInitiationEnergyFromFile("data/hairpin_initiation.data", hairpin_init);
  parsing::ParseHairpinMiscDataFromFile("data/hairpin_misc.data");
  parsing::ParseBulgeMiscDataFromFile("data/bulge_misc.data");

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              internal_1x1[a][b][c][d][e][f] = energy::MAX_E;
              for (int g = 0; g < 4; ++g) {
                internal_1x2[a][b][c][d][e][f][g] = energy::MAX_E;
                for (int h = 0; h < 4; ++h) {
                  internal_2x2[a][b][c][d][e][f][g][h] = energy::MAX_E;
                }
              }
            }
          }
        }
      }
    }
  }

  parsing::ParseInternalLoop1x1FromFile("data/internal_1x1.data");
  parsing::ParseInternalLoop1x2FromFile("data/internal_1x2.data");
  parsing::ParseInternalLoop2x2FromFile("data/internal_2x2.data");
}

base_t CharToBase(char c) {
  switch (c) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'U':
      return U;
    default:
      return -1;
  }
}

char BaseToChar(base_t b) {
  if (b < 0 || b > 4) return '?';
  return "ACGU"[b];
}

}
