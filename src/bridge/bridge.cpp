#include <constants.h>
#include "energy/energy.h"
#include "fold/fold.h"
#include "parsing.h"
#include "bridge.h"
#include "energy/structure.h"
#include "energy/energy_model.h"

#include "nn_unpaired_folder.hpp"
#include "nn_scorer.hpp"

namespace memerna {
namespace bridge {

Rnark::Rnark(const std::string& data_path) : model(data_path) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
  model.SetMLParams(93, -6, 0, 0, 999999);
}

energy_t Rnark::Efn(const secondary_t& secondary, std::string* desc) const {
  auto r = librnary::StringToPrimary(parsing::PrimaryToString(secondary.r));
  auto ss_tree = librnary::SSTree(librnary::DotBracketToMatching(parsing::PairsToDotBracket(secondary.p)));

  librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
  scorer.SetRNA(r);
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

computed_t Rnark::Fold(const primary_t& r) const {
  librnary::NNUnpairedFolder folder(model);
  folder.SetMaxTwoLoop(constants::TWOLOOP_MAX_SZ);
  folder.SetLonelyPairs(true);
  folder.SetStacking(true);
  auto librnary_primary = librnary::StringToPrimary(parsing::PrimaryToString(r));
  energy_t energy = folder.Fold(librnary_primary);
  auto pairs = parsing::DotBracketToPairs(librnary::MatchingToDotBracket(folder.Traceback()));
  return {{r, pairs}, std::vector<Ctd>(r.size(), CTD_NA), energy};
}

Rnastructure::Rnastructure(const std::string& data_path, bool use_lyngso_) :
    data(librnary::LoadDatatable(data_path)), use_lyngso(use_lyngso_) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
}

energy_t Rnastructure::Efn(const secondary_t& secondary, std::string*) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::PrimaryToString(secondary.r)),
      librnary::DotBracketToMatching(parsing::PairsToDotBracket(secondary.p))
  );
  efn2(data.get(), structure.get(), 1, true);
  return energy_t(structure->GetEnergy(1));
}

computed_t Rnastructure::Fold(const primary_t& r) const {
  return FoldAndDpTable(r, nullptr);
}

computed_t Rnastructure::FoldAndDpTable(const primary_t& r, dp_state_t* dp_state) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::PrimaryToString(r)));
  // First false here says also generate the folding itself (not just the MFE).
  // Third last parameter is whether to generate the mfe structure only -- i.e. just one.
  // Second last parameter is whether to use Lyngso or not.
  // Last parameter is for returning the DP state.
  // Add two to TWOLOOP_MAX_SZ because rnastructure bug.
  dynamic(structure.get(), data.get(), 1, 0, 0, nullptr, false, nullptr, constants::TWOLOOP_MAX_SZ + 2, true,
      !use_lyngso, dp_state);
  auto pairs = parsing::DotBracketToPairs(
      librnary::MatchingToDotBracket(librnary::StructureToMatching(*structure)));
  return {{r, pairs}, std::vector<Ctd>(r.size(), CTD_NA), energy_t(structure->GetEnergy(1))};
}

energy_t Memerna::Efn(const secondary_t& secondary, std::string* desc) const {
  computed_t computed;
  if (desc) {
    std::unique_ptr<energy::Structure> structure;
    computed = energy::ComputeEnergy(secondary, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    computed = energy::ComputeEnergy(secondary);
  }

  return computed.energy;
}

computed_t Memerna::Fold(const primary_t& r) const {
  return fold_fn(r, nullptr);
}

computed_t Memerna::FoldAndDpTable(const primary_t& r, fold::fold_state_t* fold_state) const {
  return fold_fn(r, fold_state);
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
    return std::unique_ptr<RnaPackage>(new Memerna(fold::FoldFunctionFromArgParse(argparse)));
  }
}

}
}
