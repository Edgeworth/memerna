// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_COMMON_PARSE_H_
#define BACKENDS_COMMON_PARSE_H_

#include <array>
#include <cstddef>
#include <cstring>
#include <istream>
#include <string>
#include <unordered_map>

#include "model/energy.h"
#include "util/error.h"

namespace mrna::md {

void Parse3Map(std::istream&& input, Energy (&output)[4][4][4]);

void Parse4Map(std::istream&& input, Energy (&output)[4][4][4][4]);

void Parse6Map(std::istream&& input, Energy (&output)[4][4][4][4][4][4]);

void Parse7Map(std::istream&& input, Energy (&output)[4][4][4][4][4][4][4]);

void Parse8Map(std::istream&& input, Energy (&output)[4][4][4][4][4][4][4][4]);

void ParseNMap(std::istream&& input, std::unordered_map<std::string, Energy>* output);

template <std::size_t N>
void ParseVec(std::istream&& input, Energy (&output)[N]) {
  std::memset(static_cast<void*>(&output), 0, sizeof(output));

  std::string energy;
  std::array<bool, N> seen{};
  int idx = 0;
  while (true) {
    input >> std::ws;
    if (input.peek() == std::char_traits<char>::eof()) break;
    verify(input >> idx, "expected index");
    verify(input >> energy, "expected energy for index {}", idx);
    verify(idx > 0 && idx < int(N), "out of bounds index {}", idx);
    verify(!seen[static_cast<std::size_t>(idx)], "duplicate index {}", idx);
    seen[static_cast<std::size_t>(idx)] = true;
    output[idx] = Energy::FromString(energy);
  }
  verify(input.eof(), "expected EOF");
}

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_PARSE_H_
