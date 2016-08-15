#include <constants.h>
#include "energy.h"
#include "fold.h"
#include "parsing.h"
#include "bridge.h"
#include "structure.h"

#include "nn_unpaired_folder.hpp"
#include "nn_scorer.hpp"

namespace memerna {
namespace bridge {

Rnark::Rnark(const std::string& data_path) : model(data_path) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
  model.SetMLParams(93, -6, 0, 0, 999999);
}

RnaPackage::results_t Rnark::Efn(const folded_rna_t& frna, bool verbose) const {
  auto rna = librnary::StringToPrimary(parsing::RnaToString(frna.r));
  auto ss_tree = librnary::SSTree(librnary::DotBracketToMatching(parsing::PairsToDotBracket(frna.p)));

  librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
  scorer.SetRNA(rna);
  energy_t energy;
  std::string desc;
  if (verbose) {
    auto trace = scorer.TraceExterior(ss_tree.RootSurface());
    desc = trace.Describe(' ', 0) + "\n";
    energy = trace.recursive_score;
  } else {
    energy = scorer.ScoreExterior(ss_tree.RootSurface());
  }
  return {energy, frna, desc};
}

RnaPackage::results_t Rnark::Fold(const rna_t& rna, bool verbose) const {
  librnary::NNUnpairedFolder folder(model);
  folder.SetMaxTwoLoop(constants::TWOLOOP_MAX_SZ);
  folder.SetLonelyPairs(true);
  folder.SetStacking(true);
  auto primary = librnary::StringToPrimary(parsing::RnaToString(rna));
  energy_t energy = folder.Fold(primary);
  auto pairs = parsing::DotBracketToPairs(librnary::MatchingToDotBracket(folder.Traceback()));
  return {energy, {rna, pairs}, verbose ? "verbose output not yet implemented" : ""};
}

Rnastructure::Rnastructure(const std::string& data_path, bool _use_lyngso) :
    data(librnary::LoadDatatable(data_path)), use_lyngso(_use_lyngso) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
}

RnaPackage::results_t Rnastructure::Efn(const folded_rna_t& frna, bool verbose) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::RnaToString(frna.r)),
      librnary::DotBracketToMatching(parsing::PairsToDotBracket(frna.p))
  );
  efn2(data.get(), structure.get(), 1, true);
  return {energy_t(structure->GetEnergy(1)), frna, verbose ? "verbose output not yet implemented" : ""};
}

RnaPackage::results_t Rnastructure::Fold(const rna_t& rna, bool verbose) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::RnaToString(rna)));
  // First false here says also generate the folding itself (not just the MFE).
  // Second last parameter is whether to generate the mfe structure only -- i.e. just one.
  // Last parameter is whether to use Lyngso or not.
  // Add two to TWOLOOP_MAX_SZ because rnastructure bug.
  dynamic(structure.get(), data.get(), 1, 0, 0, nullptr, false, nullptr, constants::TWOLOOP_MAX_SZ + 2, true, !use_lyngso);
  auto pairs = parsing::DotBracketToPairs(
      librnary::MatchingToDotBracket(librnary::StructureToMatching(*structure)));
  return {energy_t(structure->GetEnergy(1)), {rna, pairs}, verbose ? "verbose output not yet implemented" : ""};
}

Memerna::Memerna(const std::string& data_path) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
  LoadEnergyModelFromDataDir(data_path);
}

RnaPackage::results_t Memerna::Efn(const folded_rna_t& frna, bool verbose) const {
  energy_t energy;
  std::string desc;
  if (verbose) {
    std::unique_ptr<structure::Structure> structure;
    energy = energy::ComputeEnergy(frna, &structure);
    for (const auto& s : structure->Description()) {
      desc += s;
      desc += "\n";
    }
  } else {
    energy = energy::ComputeEnergy(frna);
  }

  return {energy, frna, desc};
}

RnaPackage::results_t Memerna::Fold(const rna_t& rna, bool verbose) const {
  energy_t energy;
  std::string desc;
  if (verbose) {
    std::unique_ptr<structure::Structure> structure;
    energy = fold::Fold(rna, &structure);
    for (const auto& s : structure->Description()) {
      desc += s;
      desc += "\n";
    }
  } else {
    energy = fold::Fold(rna);
  }
  return {energy, {rna, p}, desc};
}

}
}
