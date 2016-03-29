#include <cstdio>
#include <cassert>
#include "base.h"
#include "parsing.h"
#include "energy.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  assert(argc == 3);
  Init();
  folded_rna_t frna = parsing::ParseViennaRna(argv[1], argv[2]);
  printf("Computing energy.\n");
  printf("Computed energy: %lf\n", energy::ComputeEnergy(frna) / 10.0);
}
