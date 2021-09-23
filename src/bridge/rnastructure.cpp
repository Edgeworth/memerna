// Copyright 2016 Eliot Courtney.
#include "bridge/rnastructure.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace mrna {
namespace bridge {

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

RNAstructure::RNAstructure(const std::string& data_path, bool use_lyngso_)
    : data(LoadDatatable(data_path)), use_lyngso(use_lyngso_) {
  verify(data_path.size() && data_path.back() == '/', "invalid data path");
}

energy_t RNAstructure::Efn(const secondary_t& secondary, std::string* /*desc*/) const {
  const auto structure = LoadStructure(secondary);
  efn2(data.get(), structure.get(), 1, true);
  return energy_t(structure->GetEnergy(1));
}

computed_t RNAstructure::Fold(const primary_t& r) const {
  dp_state_t state;
  return FoldAndDpTable(r, &state);
}

computed_t RNAstructure::FoldAndDpTable(const primary_t& r, dp_state_t* dp_state) const {
  const auto structure = LoadStructure(r);
  // First false here says also generate the folding itself (not just the MFE).
  // Third last parameter is whether to generate the mfe structure only -- i.e. just one.
  // Second last parameter is whether to use Lyngso or not.
  // Last parameter is for returning the DP state.
  // Add two to TWOLOOP_MAX_SZ because rnastructure bug.
  fprintf(stderr, "DATA: %hu\n", data->auend);
  dynamic(structure.get(), data.get(), 1, 0, 0, nullptr, false, nullptr, TWOLOOP_MAX_SZ + 2, true,
      !use_lyngso, dp_state);
  return {{r, StructureToPairs(*structure)}, std::vector<Ctd>(r.size(), CTD_NA),
      energy_t(structure->GetEnergy(1))};
}

int RNAstructure::Suboptimal(
    fold::SuboptimalCallback fn, const primary_t& r, energy_t energy_delta) const {
  auto computeds = SuboptimalIntoVector(r, energy_delta);
  for (const auto& computed : computeds) fn(computed);
  return static_cast<int>(computeds.size());
}

std::vector<computed_t> RNAstructure::SuboptimalIntoVector(
    const primary_t& r, energy_t energy_delta) const {
  const auto structure = LoadStructure(r);
  // Arguments: structure, data tables, percentage delta, absolute delta, nullptr, nullptr, false
  verify(int16_t(energy_delta) == energy_delta, "energy_delta too big");
  alltrace(structure.get(), data.get(), 100, int16_t(energy_delta), nullptr, nullptr, false);
  auto p_list = StructureToMultiplePairs(*structure);
  std::vector<computed_t> computeds;
  for (int i = 0; i < static_cast<int>(p_list.size()); ++i)
    computeds.emplace_back(secondary_t(r, std::move(p_list[i])), std::vector<Ctd>(r.size(), CTD_NA),
        energy_t(structure->GetEnergy(i + 1)));
  return computeds;
}

std::pair<partition::partition_t, partition::probabilities_t> RNAstructure::Partition(
    const primary_t& r) const {
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
  const auto pfdata = std::make_unique<pfdatatable>(data.get(), scaling, T);
  const auto fce = std::make_unique<forceclass>(length);
  const auto lfce = std::make_unique<bool[]>(std::size_t(2 * length + 1));
  const auto mod = std::make_unique<bool[]>(std::size_t(2 * length + 1));

  calculatepfunction(structure.get(), pfdata.get(), nullptr, nullptr, false, nullptr, &w, &v, &wmb,
      &wl, &wlc, &wmbl, &wcoax, fce.get(), w5.get(), w3.get(), mod.get(), lfce.get());

  partition::partition_t partition = {array3d_t<penergy_t, 1>(std::size_t(length)), 0};
  partition.q = penergy_t(w5[length]);
  for (int i = 1; i <= length; ++i) {
    for (int j = i; j < length + i; ++j) {
      int adjusted = j > length ? j - length - 1 : j - 1;
      partition.p[i - 1][adjusted][0] = penergy_t(v.f(i, j));
    }
  }

  partition::probabilities_t probability(std::size_t(length), 0);
  for (int i = 0; i < length; ++i) {
    for (int j = i; j < length; ++j) {
      probability[i][j][0] = penergy_t(calculateprobability(i + 1, j + 1, &v, w5.get(),
          structure.get(), pfdata.get(), lfce.get(), mod.get(), scaling, fce.get()));
    }
  }
  return {std::move(partition), std::move(probability)};
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const primary_t& r) const {
  auto struc = std::make_unique<structure>();
  struc->SetSequence(parsing::PrimaryToString(r));
  struc->SetThermodynamicDataTable(data.get());
  return struc;
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const secondary_t& s) const {
  auto struc = LoadStructure(s.r);
  struc->AddStructure();
  for (int i = 0; i < static_cast<int>(s.p.size()); ++i) {
    if (i < s.p[i]) struc->SetPair(i + 1, s.p[i] + 1);
  }
  return struc;
}

}  // namespace bridge
}  // namespace mrna
