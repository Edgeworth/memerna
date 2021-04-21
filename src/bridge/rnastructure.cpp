// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#include "bridge/rnastructure.h"

#include "miles_rnastructure/include/alltrace.h"

namespace memerna {
namespace bridge {

namespace {

const char* FILENAMES[] = {"loop.dat", "stack.dat", "tstackh.dat", "tstacki.dat", "tloop.dat",
    "miscloop.dat", "dangle.dat", "int22.dat", "int21.dat", "coaxial.dat", "tstackcoax.dat",
    "coaxstack.dat", "tstack.dat", "tstackm.dat", "triloop.dat", "int11.dat", "hexaloop.dat",
    "tstacki23.dat", "tstacki1n.dat"};
constexpr int NUM_FILENAMES = sizeof(FILENAMES) / sizeof(FILENAMES[0]);

// Paraphrased with permission from Max Ward.
std::unique_ptr<datatable> LoadDatatable(const std::string& path) {
  auto dt = std::make_unique<datatable>();
  std::string paths[NUM_FILENAMES];
  for (int i = 0; i < NUM_FILENAMES; ++i)
    paths[i] = path + FILENAMES[i];

  opendat(paths[0].c_str(), paths[1].c_str(), paths[2].c_str(), paths[3].c_str(), paths[4].c_str(),
      paths[5].c_str(), paths[6].c_str(), paths[7].c_str(), paths[8].c_str(), paths[9].c_str(),
      paths[10].c_str(), paths[11].c_str(), paths[12].c_str(), paths[13].c_str(), paths[14].c_str(),
      paths[15].c_str(), paths[16].c_str(), paths[17].c_str(), paths[18].c_str(), dt.get());
  return dt;
}

std::unique_ptr<structure> LoadStructure(const primary_t& r) {
  auto struc = std::make_unique<structure>();
  struc->allocate(int(r.size()));
  for (int i = 0; i < int(r.size()); ++i) {
    struc->numseq[i + 1] = short(r[i] + 1);
    struc->nucs[i + 1] = BaseToChar(r[i]);
    struc->hnumber[i + 1] = short(i + 1);
  }
  return struc;
}

std::unique_ptr<structure> LoadStructure(const secondary_t& s) {
  auto struc = LoadStructure(s.r);
  struc->AddStructure();
  for (int i = 0; i < int(s.p.size()); ++i) {
    if (i < s.p[i]) struc->SetPair(i + 1, s.p[i] + 1);
  }
  return struc;
}

std::vector<int> StructureToPairs(structure& struc, int struc_num = 1) {
  std::vector<int> p(struc.GetSequenceLength(), -1);
  for (int i = 0; i < int(p.size()); ++i)
    p[i] = struc.GetPair(i + 1, struc_num) - 1;
  return p;
}

std::vector<std::vector<int>> StructureToMultiplePairs(structure& struc) {
  std::vector<std::vector<int>> ps;
  for (int i = 0; i < struc.GetNumberofStructures(); ++i)
    ps.push_back(StructureToPairs(struc, i + 1));
  return ps;
}
}

Rnastructure::Rnastructure(const std::string& data_path, bool use_lyngso_)
    : data(LoadDatatable(data_path)), use_lyngso(use_lyngso_) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
}

energy_t Rnastructure::Efn(const secondary_t& secondary, std::string* /*desc*/) const {
  const auto structure = LoadStructure(secondary);
  efn2(data.get(), structure.get(), 1, true);
  return energy_t(structure->GetEnergy(1));
}

computed_t Rnastructure::Fold(const primary_t& r) const {
  dp_state_t state;
  return FoldAndDpTable(r, state);
}

computed_t Rnastructure::FoldAndDpTable(const primary_t& r, dp_state_t& dp_state) const {
  const auto structure = LoadStructure(r);
  // First false here says also generate the folding itself (not just the MFE).
  // Third last parameter is whether to generate the mfe structure only -- i.e. just one.
  // Second last parameter is whether to use Lyngso or not.
  // Last parameter is for returning the DP state.
  // Add two to TWOLOOP_MAX_SZ because rnastructure bug.
  dynamic(structure.get(), data.get(), 1, 0, 0, nullptr, false, nullptr, TWOLOOP_MAX_SZ + 2, true,
      !use_lyngso, &dp_state);
  return {{r, StructureToPairs(*structure)}, std::vector<Ctd>(r.size(), CTD_NA),
      energy_t(structure->GetEnergy(1))};
}

int Rnastructure::Suboptimal(fold::SuboptimalCallback fn,
    const primary_t& r, energy_t energy_delta) const {
  auto computeds = SuboptimalIntoVector(r, energy_delta);
  for (const auto& computed : computeds)
    fn(computed);
  return int(computeds.size());
}

std::vector<computed_t> Rnastructure::SuboptimalIntoVector(
    const primary_t& r, energy_t energy_delta) const {
  const auto structure = LoadStructure(r);
  // Arguments: structure, data tables, percentage delta, absolute delta, nullptr, nullptr, false
  verify_expr(short(energy_delta) == energy_delta, "energy_delta too big");
  alltrace(structure.get(), data.get(), 100, short(energy_delta), nullptr, nullptr, false);
  auto p_list = StructureToMultiplePairs(*structure);
  std::vector<computed_t> computeds;
  for (int i = 0; i < int(p_list.size()); ++i)
    computeds.emplace_back(
        secondary_t(r, std::move(p_list[i])), std::vector<Ctd>(r.size(), CTD_NA),
        energy_t(structure->GetEnergy(i + 1)));
  return computeds;
}

std::pair<partition::partition_t, partition::probabilities_t> Rnastructure::Partition(
    const primary_t& r) const {
  const auto structure = LoadStructure(r);
  const int length = int(r.size());
  const PFPRECISION scaling = 1.0;  // TODO return scaling to 0.6.
  DynProgArray<PFPRECISION> w(length);
  DynProgArray<PFPRECISION> v(length);
  DynProgArray<PFPRECISION> wmb(length);
  DynProgArray<PFPRECISION> wl(length);
  DynProgArray<PFPRECISION> wmbl(length);
  DynProgArray<PFPRECISION> wcoax(length);
  const auto w5 = std::make_unique<PFPRECISION[]>(std::size_t(length + 1));
  const auto w3 = std::make_unique<PFPRECISION[]>(std::size_t(length + 2));
  const auto pfdata = std::make_unique<pfdatatable>(data.get(), scaling, T);
  const auto fce = std::make_unique<forceclass>(length);
  const auto lfce = std::make_unique<bool[]>(std::size_t(2 * length + 1));
  const auto mod = std::make_unique<bool[]>(std::size_t(2 * length + 1));

  calculatepfunction(structure.get(), pfdata.get(), nullptr, nullptr, false, nullptr,
      &w, &v, &wmb, &wl, &wmbl, &wcoax, fce.get(), w5.get(), w3.get(), mod.get(), lfce.get());

  partition::partition_t partition = {{std::size_t(length)}, 0};
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
      probability[i][j][0] = penergy_t(calculateprobability(
          i + 1, j + 1, &v, w5.get(), structure.get(),
          pfdata.get(), lfce.get(), mod.get(), scaling, fce.get()));
    }
  }
  return {std::move(partition), std::move(probability)};
}

}
}
