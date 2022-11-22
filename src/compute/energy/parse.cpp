#include "compute/energy/parse.h"

namespace mrna::energy {

void Parse3MapFromFile(const std::string& filename, Energy (&output)[4][4][4]) {
  std::ifstream f(filename);
  verify(f, "could not open file");

  std::string energy;
  while (true) {
    const auto a = CharToBase(static_cast<char>(f.get()));
    const auto b = CharToBase(static_cast<char>(f.get()));
    const auto c = CharToBase(static_cast<char>(f.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value(), "expected base");
    verify(f >> energy, "expected energy");
    output[*a][*b][*c] = Energy::FromString(energy);
  }
  verify(f.eof(), "expected EOF");
}

void Parse4MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4]) {
  std::ifstream f(filename);
  verify(f, "could not open file");

  std::string energy;
  while (true) {
    const auto a = CharToBase(static_cast<char>(f.get()));
    const auto b = CharToBase(static_cast<char>(f.get()));
    const auto c = CharToBase(static_cast<char>(f.get()));
    const auto d = CharToBase(static_cast<char>(f.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value(), "expected base");
    verify(f >> energy, "expected energy");
    output[*a][*b][*c][*d] = Energy::FromString(energy);
  }
  verify(f.eof(), "expected EOF");
}

void Parse6MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4]) {
  std::ifstream f(filename);
  verify(f, "could not open file");

  std::string energy;
  while (true) {
    const auto a = CharToBase(static_cast<char>(f.get()));
    const auto b = CharToBase(static_cast<char>(f.get()));
    const auto c = CharToBase(static_cast<char>(f.get()));
    const auto d = CharToBase(static_cast<char>(f.get()));
    const auto e = CharToBase(static_cast<char>(f.get()));
    const auto f = CharToBase(static_cast<char>(f.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value() && e.has_value() &&
            f.has_value(),
        "expected base");
    verify(f >> energy, "expected energy");
    output[*a][*b][*c][*d][*e][*f] = Energy::FromString(energy);
  }
  verify(f.eof(), "expected EOF");
}

void Parse7MapFromFile(const std::string& filename, Energy (&output)[4][4][4][4][4][4][4]) {
  std::ifstream f(filename);
  verify(f, "could not open file");

  std::string energy;
  while (true) {
    const auto a = CharToBase(static_cast<char>(f.get()));
    const auto b = CharToBase(static_cast<char>(f.get()));
    const auto c = CharToBase(static_cast<char>(f.get()));
    const auto d = CharToBase(static_cast<char>(f.get()));
    const auto e = CharToBase(static_cast<char>(f.get()));
    const auto f = CharToBase(static_cast<char>(f.get()));
    const auto g = CharToBase(static_cast<char>(f.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value() && e.has_value() &&
            f.has_value() && g.has_value(),
        "expected base");
    verify(f >> energy, "expected energy");
    output[*a][*b][*c][*d][*e][*f][*g] = Energy::FromString(energy);
  }
  verify(f.eof(), "expected EOF");
}

void Parse8MapFromFile(
    const std::string& filename, Energy (&output)[4][4][4][4][4][4][4][4]) {
  std::ifstream f(filename);
  verify(f, "could not open file");

  std::string energy;
  while (true) {
    const auto a = CharToBase(static_cast<char>(f.get()));
    const auto b = CharToBase(static_cast<char>(f.get()));
    const auto c = CharToBase(static_cast<char>(f.get()));
    const auto d = CharToBase(static_cast<char>(f.get()));
    const auto e = CharToBase(static_cast<char>(f.get()));
    const auto f = CharToBase(static_cast<char>(f.get()));
    const auto g = CharToBase(static_cast<char>(f.get()));
    const auto h = CharToBase(static_cast<char>(f.get()));
    if (!a.has_value()) break;
    verify(a.has_value() && b.has_value() && c.has_value() && d.has_value() && e.has_value() &&
            f.has_value() && g.has_value() && h.has_value(),
        "expected base");
    verify(f >> energy, "expected energy");
    output[*a][*b][*c][*d][*e][*f][*g][*h] = Energy::FromString(energy);
  }
  verify(f.eof(), "expected EOF");
}

void ParseNMapFromFile(
    const std::string& filename, std::unordered_map<std::string, Energy>* output) {
  std::ifstream f(filename);
  verify(f, "could not open file");

  std::string name;
  std::string energy;
  while (f >> name >> energy) (*output)[name] = Energy::FromString(energy);
  verify(f.eof(), "expected EOF");
}

}  // namespace mrna::energy
