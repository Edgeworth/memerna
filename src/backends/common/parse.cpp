// Copyright 2022 Eliot Courtney.
#include "backends/common/parse.h"

#include <array>
#include <cstddef>
#include <cstring>
#include <istream>
#include <string>
#include <unordered_map>

#include "model/base.h"

namespace mrna::md {

namespace {

template <std::size_t N, typename Output, typename SetEnergy>
void ParseMap(std::istream& input, Output& output, SetEnergy&& set_energy) {
  std::memset(static_cast<void*>(&output), 0, sizeof(output));

  constexpr std::size_t key_count = std::size_t{1} << (2 * N);
  std::array<bool, key_count> seen{};
  std::string key;
  std::string energy;
  while (true) {
    input >> std::ws;
    if (input.peek() == std::char_traits<char>::eof()) break;
    verify(input >> key, "expected bases");
    verify(key.size() == N, "expected {} bases", N);
    verify(input >> energy, "expected energy for {}", key);

    std::array<Base, N> bases{};
    std::size_t linear_index = 0;
    for (std::size_t i = 0; i < N; ++i) {
      const auto base = CharToBase(key[i]);
      verify(base.has_value(), "expected base");
      bases[i] = *base;
      linear_index = linear_index * 4 + *base;
    }
    verify(!seen[linear_index], "duplicate entry {}", key);
    seen[linear_index] = true;
    set_energy(bases, Energy::FromString(energy));
  }
  verify(input.eof(), "expected EOF");
}

}  // namespace

void Parse3Map(std::istream&& input, Energy (&output)[4][4][4]) {
  ParseMap<3>(
      input, output, [&](const auto& b, Energy energy) { output[b[0]][b[1]][b[2]] = energy; });
}

void Parse4Map(std::istream&& input, Energy (&output)[4][4][4][4]) {
  ParseMap<4>(input, output,
      [&](const auto& b, Energy energy) { output[b[0]][b[1]][b[2]][b[3]] = energy; });
}

void Parse6Map(std::istream&& input, Energy (&output)[4][4][4][4][4][4]) {
  ParseMap<6>(input, output,
      [&](const auto& b, Energy energy) { output[b[0]][b[1]][b[2]][b[3]][b[4]][b[5]] = energy; });
}

void Parse7Map(std::istream&& input, Energy (&output)[4][4][4][4][4][4][4]) {
  ParseMap<7>(input, output, [&](const auto& b, Energy energy) {
    output[b[0]][b[1]][b[2]][b[3]][b[4]][b[5]][b[6]] = energy;
  });
}

void Parse8Map(std::istream&& input, Energy (&output)[4][4][4][4][4][4][4][4]) {
  ParseMap<8>(input, output, [&](const auto& b, Energy energy) {
    output[b[0]][b[1]][b[2]][b[3]][b[4]][b[5]][b[6]][b[7]] = energy;
  });
}

void ParseNMap(std::istream&& input, std::unordered_map<std::string, Energy>* output) {
  output->clear();

  std::string name;
  std::string energy;
  while (true) {
    input >> std::ws;
    if (input.peek() == std::char_traits<char>::eof()) break;
    verify(input >> name, "expected name");
    verify(input >> energy, "expected energy for {}", name);
    const bool inserted = output->emplace(name, Energy::FromString(energy)).second;
    verify(inserted, "duplicate entry {}", name);
  }
  verify(input.eof(), "expected EOF");
}

}  // namespace mrna::md
