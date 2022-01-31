// Copyright 2016 Eliot Courtney.
#include "bridge/rnastructure.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"
#include "compute/subopt/subopt.h"
#include "util/error.h"

namespace mrna::bridge {

namespace {

std::unique_ptr<datatable> LoadDatatable(const std::string& path) {
  setDataPath(path.c_str());  // Set RNAstructure data path.
  auto dt = std::make_unique<datatable>();
  verify(dt->opendat(path.c_str(), "rna") == 1, "could not load RNAstructure data tables");
  return dt;
}

Secondary StructureToSecondary(structure& struc, int struc_num = 1) {
  Secondary s(struc.GetSequenceLength());
  for (int i = 0; i < static_cast<int>(s.size()); ++i) s[i] = struc.GetPair(i + 1, struc_num) - 1;
  return s;
}

std::vector<Secondary> StructureToSecondarys(structure& struc) {
  std::vector<Secondary> s;
  for (int i = 0; i < struc.GetNumberofStructures(); ++i)
    s.push_back(StructureToSecondary(struc, i + 1));
  return s;
}

}  // namespace

RNAstructure::RNAstructure(const std::string& data_path, bool use_lyngso)
    : data_(LoadDatatable(data_path)), use_lyngso_(use_lyngso) {
  verify(data_path.size() && data_path.back() == '/', "invalid data path");
  verify(data_->loadedTables, "BUG: data tables not loaded");
  verify(data_->loadedAlphabet, "BUG: alphabet not loaded");
}

energy::EnergyResult RNAstructure::Efn(Primary r, Secondary s, std::string* desc) const {
  const auto structure = LoadStructure(r, s);
  constexpr auto linear_multiloop = true;  // Use same efn calculation as DP.
  std::stringstream sstr;
  efn2(data_.get(), structure.get(), 1, linear_multiloop, desc ? &sstr : nullptr);
  if (desc) *desc = sstr.str();
  // TODO: convert ctds and structure?
  return energy::EnergyResult(structure->GetEnergy(1), Ctds(), nullptr);
}

ctx::FoldResult RNAstructure::Fold(Primary r) const {
  dp_state_t state;
  return FoldAndDpTable(std::move(r), &state);
}

ctx::FoldResult RNAstructure::FoldAndDpTable(Primary r, dp_state_t* dp_state) const {
  const auto structure = LoadStructure(r);
  constexpr auto num_tracebacks = 1;  // Number of structures to return. We just want one.
  constexpr auto percent_sort = 0;
  constexpr auto window = 0;
  constexpr auto progress = nullptr;
  constexpr auto energy_only = false;
  constexpr auto save_file = nullptr;
  constexpr auto max_twoloop = TWOLOOP_MAX_SZ + 2;  // BUG: Add two to TWOLOOP_MAX_SZ.
  constexpr auto mfe_structure_only = true;
  constexpr auto disable_coax = false;
  dynamic(structure.get(), data_.get(), num_tracebacks, percent_sort, window, progress, energy_only,
      save_file, max_twoloop, mfe_structure_only, !use_lyngso_, disable_coax, dp_state);
  // TODO: convert dp tables, ext, ctds?, delete this function and move all to Fold.
  return {.mfe = {.dp{}, .ext{}, .energy = Energy(structure->GetEnergy(1))},
      .tb = tb::TracebackResult(StructureToSecondary(*structure), Ctds())};
}

int RNAstructure::Suboptimal(subopt::SuboptCallback fn, Primary r, Energy delta) const {
  auto res = SuboptimalIntoVector(std::move(r), delta);
  for (const auto& r : res) fn(r);
  return static_cast<int>(res.size());
}

std::vector<subopt::SuboptResult> RNAstructure::SuboptimalIntoVector(
    Primary r, Energy delta) const {
  const auto structure = LoadStructure(r);
  // Arguments: structure, data tables, percentage delta, absolute delta, nullptr, nullptr, false
  verify(int16_t(delta) == delta, "delta too big");
  alltrace(structure.get(), data_.get(), 100, int16_t(delta), nullptr, nullptr, false);
  auto s_list = StructureToSecondarys(*structure);
  std::vector<subopt::SuboptResult> res;
  for (int i = 0; i < static_cast<int>(s_list.size()); ++i) {
    // TODO: Convert CTDs?
    res.emplace_back(subopt::SuboptResult(
        tb::TracebackResult(std::move(s_list[i]), Ctds()), Energy(structure->GetEnergy(i + 1))));
  }
  return res;
}

part::PartResult RNAstructure::Partition(Primary r) const {
  const auto structure = LoadStructure(r);
  const int length = static_cast<int>(r.size());
  const PFPRECISION scaling = 1.0;  // TODO return scaling to 0.6.
  DynProgArray<PFPRECISION> w(length);
  DynProgArray<PFPRECISION> v(length);
  DynProgArray<PFPRECISION> wmb(length);
  DynProgArray<PFPRECISION> wl(length);
  DynProgArray<PFPRECISION> wlc(length);
  DynProgArray<PFPRECISION> wmbl(length);
  DynProgArray<PFPRECISION> wcoax(length);
  const auto w5 = std::make_unique<PFPRECISION[]>(length + 1);
  const auto w3 = std::make_unique<PFPRECISION[]>(length + 2);
  const auto pfdata = std::make_unique<pfdatatable>(data_.get(), scaling, T);
  const auto fce = std::make_unique<forceclass>(length);
  const auto lfce = std::make_unique<bool[]>(2 * length + 1);
  const auto mod = std::make_unique<bool[]>(2 * length + 1);

  calculatepfunction(structure.get(), pfdata.get(), nullptr, nullptr, false, nullptr, &w, &v, &wmb,
      &wl, &wlc, &wmbl, &wcoax, fce.get(), w5.get(), w3.get(), mod.get(), lfce.get());

  part::Part partition = {BoltzSums(length, 0), 0};
  partition.q = BoltzEnergy(w5[length]);
  for (int i = 1; i <= length; ++i) {
    for (int j = i; j < length + i; ++j) {
      int adjusted = j > length ? j - length - 1 : j - 1;
      partition.p[i - 1][adjusted] = BoltzEnergy(v.f(i, j));
    }
  }

  BoltzProbs prob(length, 0);
  for (int i = 0; i < length; ++i) {
    for (int j = i; j < length; ++j) {
      prob[i][j] = BoltzEnergy(calculateprobability(i + 1, j + 1, &v, w5.get(), structure.get(),
          pfdata.get(), lfce.get(), mod.get(), scaling, fce.get()));
    }
  }
  // TODO: Convert tables?
  return {.dp{}, .ext{}, .part = std::move(partition), .prob = std::move(prob)};
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const Primary& r) const {
  auto struc = std::make_unique<structure>();
  struc->SetThermodynamicDataTable(data_.get());
  struc->SetSequence(r.ToString());
  verify(struc->GetSequenceLength() == static_cast<int>(r.size()), "BUG: structure not loaded");
  return struc;
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const Primary& r, const Secondary& s) const {
  auto struc = LoadStructure(r);
  struc->AddStructure();
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    if (i < s[i]) struc->SetPair(i + 1, s[i] + 1);
  }
  return struc;
}

}  // namespace mrna::bridge
