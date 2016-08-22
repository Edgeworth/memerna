#include <constants.h>
#include "energy/energy.h"
#include "fold/fold.h"
#include "parsing.h"
#include "bridge.h"
#include "energy/structure.h"

#include "nn_unpaired_folder.hpp"
#include "nn_scorer.hpp"

namespace memerna {
namespace bridge {

Rnark::Rnark(const std::string& data_path) : model(data_path) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
  model.SetMLParams(93, -6, 0, 0, 999999);
}

energy_t Rnark::Efn(const folded_rna_t& frna, std::string* desc) const {
  auto rna = librnary::StringToPrimary(parsing::RnaToString(frna.r));
  auto ss_tree = librnary::SSTree(librnary::DotBracketToMatching(parsing::PairsToDotBracket(frna.p)));

  librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
  scorer.SetRNA(rna);
  energy_t energy;
  if (desc) {
    auto trace = scorer.TraceExterior(ss_tree.RootSurface());
    *desc = trace.Describe(' ', 0) + "\n";
    energy = trace.recursive_score;
  } else {
    energy = scorer.ScoreExterior(ss_tree.RootSurface());
  }
  return energy;
}

folded_rna_t Rnark::Fold(const rna_t& rna) const {
  librnary::NNUnpairedFolder folder(model);
  folder.SetMaxTwoLoop(constants::TWOLOOP_MAX_SZ);
  folder.SetLonelyPairs(true);
  folder.SetStacking(true);
  auto primary = librnary::StringToPrimary(parsing::RnaToString(rna));
  energy_t energy = folder.Fold(primary);
  auto pairs = parsing::DotBracketToPairs(librnary::MatchingToDotBracket(folder.Traceback()));
  return {rna, pairs, energy};
}

Rnastructure::Rnastructure(const std::string& data_path, bool use_lyngso_) :
    data(librnary::LoadDatatable(data_path)), use_lyngso(use_lyngso_) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
}

energy_t Rnastructure::Efn(const folded_rna_t& frna, std::string*) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::RnaToString(frna.r)),
      librnary::DotBracketToMatching(parsing::PairsToDotBracket(frna.p))
  );
  efn2(data.get(), structure.get(), 1, true);
  return energy_t(structure->GetEnergy(1));
}

folded_rna_t Rnastructure::Fold(const rna_t& rna) const {
  return FoldAndDpTable(rna, nullptr);
}

folded_rna_t Rnastructure::FoldAndDpTable(const rna_t& rna, dp_state_t* dp_state) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::RnaToString(rna)));
  // First false here says also generate the folding itself (not just the MFE).
  // Third last parameter is whether to generate the mfe structure only -- i.e. just one.
  // Second last parameter is whether to use Lyngso or not.
  // Last parameter is for returning the DP state.
  // Add two to TWOLOOP_MAX_SZ because rnastructure bug.
  dynamic(structure.get(), data.get(), 1, 0, 0, nullptr, false, nullptr, constants::TWOLOOP_MAX_SZ + 2, true,
      !use_lyngso, dp_state);
  auto pairs = parsing::DotBracketToPairs(
      librnary::MatchingToDotBracket(librnary::StructureToMatching(*structure)));
  return {rna, pairs, energy_t(structure->GetEnergy(1))};
}

Memerna::Memerna(const std::string& data_path, const fold::fold_fn_t* fold_fn_) : fold_fn(fold_fn_) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
  LoadEnergyModelFromDataDir(data_path);
}

Memerna::Memerna(Memerna&& meme) {
  fold_fn = meme.fold_fn;
}

energy_t Memerna::Efn(const folded_rna_t& frna, std::string* desc) const {
  energy_t energy;
  if (desc) {
    std::unique_ptr<structure::Structure> structure;
    energy = energy::ComputeEnergy(frna, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    energy = energy::ComputeEnergy(frna);
  }

  return energy;
}

folded_rna_t Memerna::Fold(const rna_t& rna) const {
  return fold_fn(rna, nullptr);
}

folded_rna_t Memerna::FoldAndDpTable(const rna_t& rna, fold::fold_state_t* fold_state) const {
  return fold_fn(rna, fold_state);
}

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse) {
  verify_expr(
      argparse.HasFlag("r") + argparse.HasFlag("m") + argparse.HasFlag("k") == 1,
      "require exactly one package flag\n%s", argparse.Usage().c_str());
  if (argparse.HasFlag("r")) {
    return std::unique_ptr<RnaPackage>(new Rnastructure("extern/rnark/data_tables/", false));
  } else if (argparse.HasFlag("m")) {
    return std::unique_ptr<RnaPackage>(new Rnark("extern/rnark/data_tables/"));
  } else {
    return std::unique_ptr<RnaPackage>(new Memerna("data/", FoldFunctionFromArgParse(argparse)));
  }
}

}
}
