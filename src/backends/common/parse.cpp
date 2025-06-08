// Copyright 2022 Eliot Courtney.
#include "backends/common/parse.h"

#include <optional>
#include <string>
#include <unordered_map>

#include "model/base.h"

namespace mrna::md {

void Parse3MapFromFile(const std::string& filename, Energy (&output)[4][4][4]) {
  std::ifstream fp(filename);
  verify(fp, "could not open file: {}", filename);

  std::string energy;
  while (true) {
    fp >> std::ws;
    const auto a = CharToBase(static_cast<char>(fp.get()));
    const auto b = CharToBase(static_cast<char>(fp.get()));
    const auto c = CharToBase(static_cast<char>(fp.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value(), "expected base");
    verify(fp >> energy, "expected energy");
    output[*a][*b][*c] = Energy::FromString(energy);
  }
  verify(fp.eof(), "expected EOF");
}

void Parse4MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4]) {
  std::ifstream fp(filename);
  verify(fp, "could not open file: {}", filename);

  std::string energy;
  while (true) {
    fp >> std::ws;
    const auto a = CharToBase(static_cast<char>(fp.get()));
    const auto b = CharToBase(static_cast<char>(fp.get()));
    const auto c = CharToBase(static_cast<char>(fp.get()));
    const auto d = CharToBase(static_cast<char>(fp.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value(), "expected base");
    verify(fp >> energy, "expected energy");
    output[*a][*b][*c][*d] = Energy::FromString(energy);
  }
  verify(fp.eof(), "expected EOF");
}

void Parse6MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4]) {
  std::ifstream fp(filename);
  verify(fp, "could not open file: {}", filename);

  std::string energy;
  while (true) {
    fp >> std::ws;
    const auto a = CharToBase(static_cast<char>(fp.get()));
    const auto b = CharToBase(static_cast<char>(fp.get()));
    const auto c = CharToBase(static_cast<char>(fp.get()));
    const auto d = CharToBase(static_cast<char>(fp.get()));
    const auto e = CharToBase(static_cast<char>(fp.get()));
    const auto f = CharToBase(static_cast<char>(fp.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value() && e.has_value() &&
            f.has_value(),
        "expected base");
    verify(fp >> energy, "expected energy");
    output[*a][*b][*c][*d][*e][*f] = Energy::FromString(energy);
  }
  verify(fp.eof(), "expected EOF");
}

void Parse7MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4][4]) {
  std::ifstream fp(filename);
  verify(fp, "could not open file: {}", filename);

  std::string energy;
  while (true) {
    fp >> std::ws;
    const auto a = CharToBase(static_cast<char>(fp.get()));
    const auto b = CharToBase(static_cast<char>(fp.get()));
    const auto c = CharToBase(static_cast<char>(fp.get()));
    const auto d = CharToBase(static_cast<char>(fp.get()));
    const auto e = CharToBase(static_cast<char>(fp.get()));
    const auto f = CharToBase(static_cast<char>(fp.get()));
    const auto g = CharToBase(static_cast<char>(fp.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value() && e.has_value() &&
            f.has_value() && g.has_value(),
        "expected base");
    verify(fp >> energy, "expected energy");
    output[*a][*b][*c][*d][*e][*f][*g] = Energy::FromString(energy);
  }
  verify(fp.eof(), "expected EOF");
}

void Parse8MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4][4][4]) {
  std::ifstream fp(filename);
  verify(fp, "could not open file: {}", filename);

  std::string energy;
  while (true) {
    fp >> std::ws;
    const auto a = CharToBase(static_cast<char>(fp.get()));
    const auto b = CharToBase(static_cast<char>(fp.get()));
    const auto c = CharToBase(static_cast<char>(fp.get()));
    const auto d = CharToBase(static_cast<char>(fp.get()));
    const auto e = CharToBase(static_cast<char>(fp.get()));
    const auto f = CharToBase(static_cast<char>(fp.get()));
    const auto g = CharToBase(static_cast<char>(fp.get()));
    const auto h = CharToBase(static_cast<char>(fp.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value() && e.has_value() &&
            f.has_value() && g.has_value() && h.has_value(),
        "expected base");
    verify(fp >> energy, "expected energy");
    output[*a][*b][*c][*d][*e][*f][*g][*h] = Energy::FromString(energy);
  }
  verify(fp.eof(), "expected EOF");
}

void ParseNMapFromFile(
    const std::string& filename, std::unordered_map<std::string, Energy>* output) {
  std::ifstream fp(filename);
  verify(fp, "could not open file: {}", filename);

  std::string name;
  std::string energy;
  while (fp >> name >> energy) (*output)[name] = Energy::FromString(energy);
  verify(fp.eof(), "expected EOF");
}

}  // namespace mrna::md
