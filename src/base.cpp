#include "base.h"
#include "parsing.h"

namespace memerna {

void Init() {
  // Stacking interaction data.
  parsing::Parse2x2FromFile("data/stacking.data", stacking_e);

  // Terminal mismatch data.
  parsing::Parse2x2FromFile("data/terminal.data", terminal_e);

  // Hairpin data.
  parsing::ParseMapFromFile("data/hairpin.data", hairpin_e);
  parsing::ParseInitiationEnergyFromFile("data/hairpin_initiation.data", hairpin_init);
  parsing::ParseHairpinMiscDataFromFile("data/hairpin_misc.data");

  // Bulge loop data.
  parsing::ParseInitiationEnergyFromFile("data/bulge_initiation.data", bulge_init);
  parsing::ParseBulgeMiscDataFromFile("data/bulge_misc.data");

  // Internal loop data.
  parsing::ParseInitiationEnergyFromFile("data/internal_initiation.data", internal_init);
  parsing::ParseInternalLoopMiscDataFromFile("data/internal_misc.data");
  parsing::ParseInternalLoop1x1FromFile("data/internal_1x1.data");
  parsing::ParseInternalLoop1x2FromFile("data/internal_1x2.data");
  parsing::ParseInternalLoop2x2FromFile("data/internal_2x2.data");
  parsing::Parse2x2FromFile("data/internal_2x3_mismatch.data", internal_2x3_mismatch);
  parsing::Parse2x2FromFile("data/internal_other_mismatch.data", internal_other_mismatch);

  // Multiloop data.
  parsing::ParseMultiloopMiscDataFromFile("data/multiloop_misc.data");
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
