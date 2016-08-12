#include <cstdio>
#include <iostream>
#include "base.h"
#include "parsing.h"
#include "energy.h"

using namespace memerna;

int main() {
  LoadEnergyModelFromDataDir();

  int idx = 0;
  while (1) {
    std::string name, seq, db;
    getline(std::cin, name);
    getline(std::cin, seq);
    getline(std::cin, db);
    if (!std::cin) break;
    auto frna = parsing::ParseDotBracketRna(seq, db);
    printf("RNA %d\nEnergy: %d\n", idx++, energy::ComputeEnergy(frna));
  }
}
