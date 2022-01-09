// Copyright 2016 E.
#include "bridge/rnastructure.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "util/macros.h"

namespace mrna::bridge {

namespace {

std::unique_ptr<datatable> LoadDatatable(const std::string& path) {
  setDataPath(path.c_str());  // Set RNAstructure data path.
  auto dt = std::make_unique<datatable>();
  verify(dt->opendat(path.c_str(), "rna") == 1, "could not load RNAstructure data tables");
  return dt;
}

std::vector<int> StructureToPairs(structure& struc, int struc_num = 1) {
  std::vector<int> p(struc.GetSequenceLength(), -1);
  for (int i = 0; i < static_cast<int>(p.size()); ++i) p[i] = struc.GetPair(i + 1, struc_num) - 1;
  return p;
}

std::vector<std::vector<int>> StructureToMultiplePairs(structure& struc) {
  std::vector<std::vector<int>> ps;
  for (int i = 0; i < struc.GetNumberofStructures(); ++i)
    ps.push_back(StructureToPairs(struc, i + 1));
  return ps;
}

}  // namespace

RNAstructure::RNAstructure(const std::string& data_path, bool use_lyngso)
    : data_(LoadDatatable(data_path)), use_lyngso_(use_lyngso) {
  verify(data_path.size() && data_path.back() == '/', "invalid data path");
  verify(data_->loadedTables, "BUG: data tables not loaded");
  verify(data_->loadedAlphabet, "BUG: alphabet not loaded");
}

Energy RNAstructure::Efn(const Secondary& secondary, std::string* desc) const {
  const auto structure = LoadStructure(secondary);
  constexpr auto linear_multiloop = true;  // Use same efn calculation as DP.
  std::stringstream sstr;
  efn2(data_.get(), structure.get(), 1, linear_multiloop, desc ? &sstr : nullptr);
  if (desc) *desc = sstr.str();
  return Energy(structure->GetEnergy(1));
}

Computed RNAstructure::Fold(const Primary& r) const {
  dp_state_t state;
  return FoldAndDpTable(r, &state);
}

Computed RNAstructure::FoldAndDpTable(const Primary& r, dp_state_t* dp_state) const {
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
  return {{r, StructureToPairs(*structure)}, std::vector<Ctd>(r.size(), CTD_NA),
      Energy(structure->GetEnergy(1))};
}

int RNAstructure::Suboptimal(
    subopt::SuboptimalCallback fn, const Primary& r, Energy energy_delta) const {
  auto computeds = SuboptimalIntoVector(r, energy_delta);
  for (const auto& computed : computeds) fn(computed);
  return static_cast<int>(computeds.size());
}

std::vector<Computed> RNAstructure::SuboptimalIntoVector(
    const Primary& r, Energy energy_delta) const {
  const auto structure = LoadStructure(r);
  // Arguments: structure, data tables, percentage delta, absolute delta, nullptr, nullptr, false
  verify(int16_t(energy_delta) == energy_delta, "energy_delta too big");
  alltrace(structure.get(), data_.get(), 100, int16_t(energy_delta), nullptr, nullptr, false);
  auto p_list = StructureToMultiplePairs(*structure);
  std::vector<Computed> computeds;
  for (int i = 0; i < static_cast<int>(p_list.size()); ++i)
    computeds.emplace_back(Secondary(r, std::move(p_list[i])), std::vector<Ctd>(r.size(), CTD_NA),
        Energy(structure->GetEnergy(i + 1)));
  return computeds;
}

std::pair<partition::Partition, Probabilities> RNAstructure::Partition(const Primary& r) const {
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
  const auto w5 = std::make_unique<PFPRECISION[]>(std::size_t(length + 1));
  const auto w3 = std::make_unique<PFPRECISION[]>(std::size_t(length + 2));
  const auto pfdata = std::make_unique<pfdatatable>(data_.get(), scaling, T);
  const auto fce = std::make_unique<forceclass>(length);
  const auto lfce = std::make_unique<bool[]>(std::size_t(2 * length + 1));
  const auto mod = std::make_unique<bool[]>(std::size_t(2 * length + 1));

  calculatepfunction(structure.get(), pfdata.get(), nullptr, nullptr, false, nullptr, &w, &v, &wmb,
      &wl, &wlc, &wmbl, &wcoax, fce.get(), w5.get(), w3.get(), mod.get(), lfce.get());

  partition::Partition partition = {Array3D<BoltzEnergy, 1>(std::size_t(length)), 0};
  partition.q = BoltzEnergy(w5[length]);
  for (int i = 1; i <= length; ++i) {
    for (int j = i; j < length + i; ++j) {
      int adjusted = j > length ? j - length - 1 : j - 1;
      partition.p[i - 1][adjusted][0] = BoltzEnergy(v.f(i, j));
    }
  }

  Probabilities probability(std::size_t(length), 0);
  for (int i = 0; i < length; ++i) {
    for (int j = i; j < length; ++j) {
      probability[i][j][0] = BoltzEnergy(calculateprobability(i + 1, j + 1, &v, w5.get(),
          structure.get(), pfdata.get(), lfce.get(), mod.get(), scaling, fce.get()));
    }
  }
  return {std::move(partition), std::move(probability)};
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const Primary& r) const {
  auto struc = std::make_unique<structure>();
  struc->SetThermodynamicDataTable(data_.get());
  struc->SetSequence(PrimaryToString(r));
  verify(struc->GetSequenceLength() == static_cast<int>(r.size()), "BUG: structure not loaded");
  return struc;
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const Secondary& s) const {
  auto struc = LoadStructure(s.r);
  struc->AddStructure();
  for (int i = 0; i < static_cast<int>(s.p.size()); ++i) {
    if (i < s.p[i]) struc->SetPair(i + 1, s.p[i] + 1);
  }
  return struc;
}

}  // namespace mrna::bridge
