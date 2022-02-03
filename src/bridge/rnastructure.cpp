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
#include "miles_rnastructure/include/stochastic.h"
#include "util/error.h"

namespace mrna::bridge {

namespace {

std::unique_ptr<datatable> LoadDatatable(const std::string& path) {
  setDataPath(path.c_str());  // Set RNAstructure data path.
  auto dt = std::make_unique<datatable>();
  verify(dt->opendat(path.c_str(), "rna") == 1, "could not load RNAstructure data tables");
  return dt;
}

Secondary StructureToSecondary(const structure& struc, int struc_num = 1) {
  Secondary s(struc.GetSequenceLength());
  for (int i = 0; i < static_cast<int>(s.size()); ++i) s[i] = struc.GetPair(i + 1, struc_num) - 1;
  return s;
}

std::vector<Secondary> StructureToSecondarys(const structure& struc) {
  std::vector<Secondary> s;
  for (int i = 0; i < struc.GetNumberofStructures(); ++i)
    s.push_back(StructureToSecondary(struc, i + 1));
  return s;
}

std::vector<subopt::SuboptResult> StructureToSuboptVector(const structure& struc) {
  auto s_list = StructureToSecondarys(struc);
  std::vector<subopt::SuboptResult> res;
  for (int i = 0; i < static_cast<int>(s_list.size()); ++i) {
    // TODO: Convert CTDs?
    res.emplace_back(subopt::SuboptResult(
        Energy(struc.GetEnergy(i + 1)), tb::TracebackResult(std::move(s_list[i]), Ctds())));
  }
  return res;
}

struct PartitionState {
  const PFPRECISION scaling = 1.0;  // TODO return scaling to 0.6.
  std::unique_ptr<DynProgArray<PFPRECISION>> w;
  std::unique_ptr<DynProgArray<PFPRECISION>> v;
  std::unique_ptr<DynProgArray<PFPRECISION>> wmb;
  std::unique_ptr<DynProgArray<PFPRECISION>> wl;
  std::unique_ptr<DynProgArray<PFPRECISION>> wlc;
  std::unique_ptr<DynProgArray<PFPRECISION>> wmbl;
  std::unique_ptr<DynProgArray<PFPRECISION>> wcoax;
  std::unique_ptr<PFPRECISION[]> w5;
  std::unique_ptr<PFPRECISION[]> w3;
  std::unique_ptr<pfdatatable> pfdata;
  std::unique_ptr<forceclass> fce;
  std::unique_ptr<bool[]> lfce;
  std::unique_ptr<bool[]> mod;
  PartitionState(int N, datatable* data)
      : w(new DynProgArray<PFPRECISION>(N)), v(new DynProgArray<PFPRECISION>(N)),
        wmb(new DynProgArray<PFPRECISION>(N)), wl(new DynProgArray<PFPRECISION>(N)),
        wlc(new DynProgArray<PFPRECISION>(N)), wmbl(new DynProgArray<PFPRECISION>(N)),
        wcoax(new DynProgArray<PFPRECISION>(N)), w5(new PFPRECISION[N + 1]),
        w3(new PFPRECISION[N + 2]), pfdata(new pfdatatable(data, scaling, T)),
        fce(new forceclass(N)), lfce(new bool[2 * N + 1]), mod(new bool[2 * N + 1]) {}
};

PartitionState RunPartition(structure* struc, datatable* data) {
  PartitionState state(struc->GetSequenceLength(), data);

  calculatepfunction(struc, state.pfdata.get(), nullptr, nullptr, false, nullptr, state.w.get(),
      state.v.get(), state.wmb.get(), state.wl.get(), state.wlc.get(), state.wmbl.get(),
      state.wcoax.get(), state.fce.get(), state.w5.get(), state.w3.get(), state.mod.get(),
      state.lfce.get());
  return state;
}

}  // namespace

RNAstructure::RNAstructure(const std::string& data_path, bool use_lyngso)
    : data_(LoadDatatable(data_path)), use_lyngso_(use_lyngso) {
  verify(data_path.size() && data_path.back() == '/', "invalid data path");
  verify(data_->loadedTables, "BUG: data tables not loaded");
  verify(data_->loadedAlphabet, "BUG: alphabet not loaded");
}

energy::EnergyResult RNAstructure::Efn(
    const Primary& r, const Secondary& s, std::string* desc) const {
  const auto structure = LoadStructure(r, s);
  constexpr auto linear_multiloop = true;  // Use same efn calculation as DP.
  std::stringstream sstr;
  efn2(data_.get(), structure.get(), 1, linear_multiloop, desc ? &sstr : nullptr);
  if (desc) *desc = sstr.str();
  // TODO: convert ctds and structure?
  return energy::EnergyResult(structure->GetEnergy(1), Ctds(), nullptr);
}

ctx::FoldResult RNAstructure::Fold(const Primary& r) const {
  dp_state_t state;
  return FoldAndDpTable(r, &state);
}

ctx::FoldResult RNAstructure::FoldAndDpTable(const Primary& r, dp_state_t* dp_state) const {
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

int RNAstructure::Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy delta) const {
  auto res = SuboptimalIntoVector(r, delta);
  for (const auto& subopt : res) fn(subopt);
  return static_cast<int>(res.size());
}

std::vector<subopt::SuboptResult> RNAstructure::SuboptimalIntoVector(
    const Primary& r, Energy delta) const {
  const auto structure = LoadStructure(r);
  // Arguments: structure, data tables, percentage delta, absolute delta, nullptr, nullptr, false
  verify(int16_t(delta) == delta, "delta too big");
  alltrace(structure.get(), data_.get(), 100, int16_t(delta), nullptr, nullptr, false);
  return StructureToSuboptVector(*structure);
}

part::PartResult RNAstructure::Partition(const Primary& r) const {
  const auto structure = LoadStructure(r);
  auto state = RunPartition(structure.get(), data_.get());
  const int N = static_cast<int>(r.size());

  part::Part part = {BoltzSums(N, 0), 0};
  part.q = BoltzEnergy(state.w5[N]);
  for (int i = 1; i <= N; ++i) {
    for (int j = i; j < N + i; ++j) {
      int adjusted = j > N ? j - N - 1 : j - 1;
      part.p[i - 1][adjusted] = BoltzEnergy(state.v->f(i, j));
    }
  }

  BoltzProbs prob(N, 0);
  for (int i = 0; i < N; ++i) {
    for (int j = i; j < N; ++j) {
      prob[i][j] = BoltzEnergy(calculateprobability(i + 1, j + 1, state.v.get(), state.w5.get(),
          structure.get(), state.pfdata.get(), state.lfce.get(), state.mod.get(), state.scaling,
          state.fce.get()));
    }
  }
  // TODO: Convert tables?
  return {.dp{}, .ext{}, .part{std::move(part)}, .prob{std::move(prob)}};
}

std::vector<subopt::SuboptResult> RNAstructure::StochasticSampleIntoVector(
    const Primary& r, int num_samples) const {
  const auto structure = LoadStructure(r);
  auto state = RunPartition(structure.get(), data_.get());
  stochastictraceback(state.w.get(), state.wmb.get(), state.wmbl.get(), state.wcoax.get(),
      state.wl.get(), state.wlc.get(), state.v.get(), state.fce.get(), state.w3.get(),
      state.w5.get(), state.scaling, state.lfce.get(), state.mod.get(), state.pfdata.get(),
      num_samples, structure.get());
  return StructureToSuboptVector(*structure);
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

RNAstructure RNAstructure::FromArgParse(const ArgParse& args) {
  return RNAstructure(args.Get(OPT_RNASTRUCTURE_DATA), false);
}

}  // namespace mrna::bridge
