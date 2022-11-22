#ifndef COMPUTE_ENERGY_PARSE_H_
#define COMPUTE_ENERGY_PARSE_H_

#include <fstream>
#include <unordered_map>

#include "model/energy.h"

namespace mrna::energy {

void Parse3MapFromFile(const std::string& filename, Energy (&output)[4][4][4]);

void Parse4MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4]);

void Parse6MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4]);

void Parse7MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4][4]);

void Parse8MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4][4][4]);

void ParseNMapFromFile(
    const std::string& filename, std::unordered_map<std::string, Energy>* output);

template <std::size_t N>
void ParseVecFromFile(const std::string& filename, Energy (&output)[N]) {
  std::ifstream f(filename);
  verify(f, "could not open file");

  std::string energy;
  int idx = 0;
  while (f >> idx >> energy) {
    verify(idx < N, "out of bounds index in %s", filename.c_str());
    output[idx] = Energy::FromString(energy);
  }
  verify(f.eof(), "expected EOF");
}

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_PARSE_H_
