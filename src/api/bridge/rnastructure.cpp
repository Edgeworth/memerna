// Copyright 2016 Eliot Courtney.
#include "api/bridge/rnastructure.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/trace/trace.h"
#include "rnastructure_bridge/include/alltrace.h"
#include "rnastructure_bridge/include/stochastic.h"
#include "util/error.h"

namespace mrna::bridge {

namespace {

std::unique_ptr<datatable> LoadDatatable(const std::string& path) {
  verify(!path.empty() && path.back() == '/', "invalid data path");
  setDataPath(path.c_str());  // Set RNAstructure data path.
  auto dt = std::make_unique<datatable>();
  verify(dt->opendat(path.c_str(), "rna") == 1, "could not load RNAstructure data tables");
  return dt;
}

Secondary StructureToSecondary(const structure& struc, int struc_num = 1) {
  Secondary s(struc.GetSequenceLength());
  for (int i = 0; i < s.size(); ++i) {
    int pair = struc.GetPair(i + 1, struc_num);
    s[i] = pair == 0 ? INVALID_INDEX : As<Index>(pair - 1);
  }
  return s;
}

std::vector<Secondary> StructureToSecondarys(const structure& struc) {
  std::vector<Secondary> s;
  s.reserve(struc.GetNumberofStructures());
  for (int i = 0; i < struc.GetNumberofStructures(); ++i)
    s.push_back(StructureToSecondary(struc, i + 1));
  return s;
}

std::vector<subopt::SuboptResult> StructureToSuboptVector(const structure& struc) {
  auto s_list = StructureToSecondarys(struc);
  std::vector<subopt::SuboptResult> res;
  res.reserve(static_cast<int>(s_list.size()));
  for (int i = 0; i < static_cast<int>(s_list.size()); ++i) {
    // TODO(2): Convert CTDs?
    res.emplace_back(RNAstructure::ToEnergy(struc.GetEnergy(i + 1)),
        trace::TraceResult(std::move(s_list[i]), Ctds()));
  }
  return res;
}

struct PfnState {
  const PFPRECISION scaling = 1.0;  // TODO(0) return scaling to 0.6.
  DynProgArray<PFPRECISION> w;
  DynProgArray<PFPRECISION> v;
  DynProgArray<PFPRECISION> wmb;
  DynProgArray<PFPRECISION> wl;
  DynProgArray<PFPRECISION> wlc;
  DynProgArray<PFPRECISION> wmbl;
  DynProgArray<PFPRECISION> wcoax;
  std::unique_ptr<PFPRECISION[]> w5;
  std::unique_ptr<PFPRECISION[]> w3;
  std::unique_ptr<pfdatatable> pfdata;
  std::unique_ptr<forceclass> fce;
  std::unique_ptr<bool[]> lfce;
  std::unique_ptr<bool[]> mod;

  PfnState(int N, datatable* data)
      : w(N), v(N), wmb(N), wl(N), wlc(N), wmbl(N), wcoax(N), w5(new PFPRECISION[N + 1]()),
        w3(new PFPRECISION[N + 2]()), pfdata(new pfdatatable(data, scaling, double(T))),
        fce(new forceclass(N)), lfce(new bool[2 * N + 1]()), mod(new bool[2 * N + 1]()) {}
};

PfnState RunPfn(structure* struc, datatable* data) {
  PfnState state(struc->GetSequenceLength(), data);

  calculatepfunction(struc, state.pfdata.get(), nullptr, nullptr, false, nullptr, &state.w,
      &state.v, &state.wmb, &state.wl, &state.wlc, &state.wmbl, &state.wcoax, state.fce.get(),
      state.w5.get(), state.w3.get(), state.mod.get(), state.lfce.get());
  return state;
}

void VerifySupported(
    std::optional<MfeAlg> alg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf) {
  verify(!alg.has_value(), "Bridge for RNAstructure does not support algorithm selection");
  verify(pf.Empty(), "Bridge for RNAstructure does not support pseudofree energies");
  verify(cfg.lonely_pairs == erg::EnergyCfg::LonelyPairs::HEURISTIC,
      "Bridge for RNAstructure does not support lonely_pairs != HEURISTIC");
  verify(
      cfg.ctd == erg::EnergyCfg::Ctd::ALL, "Bridge for RNAstructure does not support ctd != ALL");
}

void VerifySupported(
    std::optional<PfnAlg> alg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf) {
  verify(!alg.has_value(), "Bridge for RNAstructure does not support algorithm selection");
  verify(pf.Empty(), "Bridge for RNAstructure does not support pseudofree energies");
  verify(cfg.lonely_pairs == erg::EnergyCfg::LonelyPairs::HEURISTIC,
      "Bridge for RNAstructure does not support lonely_pairs != HEURISTIC");
  verify(
      cfg.ctd == erg::EnergyCfg::Ctd::ALL, "Bridge for RNAstructure does not support ctd != ALL");
}

void VerifySupported(const trace::TraceCfg& trace_cfg) {
  verify(!trace_cfg.random, "Bridge for RNAstructure does not support random trace");
}

void VerifySupported(const subopt::SuboptCfg& subopt_cfg) {
  verify(subopt_cfg.strucs == subopt::SuboptCfg::MAX_STRUCTURES,
      "Bridge for RNAstructure does not support strucs limit");
  verify(subopt_cfg.time_secs < 0, "Bridge for RNAstructure does not support time limit");
}

}  // namespace

RNAstructure::RNAstructure(const std::string& data_path, bool use_lyngso)
    : data_(LoadDatatable(data_path)), use_lyngso_(use_lyngso) {
  verify(data_->loadedTables, "BUG: data tables not loaded");
  verify(data_->loadedAlphabet, "BUG: alphabet not loaded");
}

erg::EnergyResult RNAstructure::Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const Ctds* given_ctd, bool /*build_structure*/) const {
  verify(given_ctd == nullptr, "Bridge for RNAstructure does not support given_ctd");
  VerifySupported(std::optional<MfeAlg>{}, cfg, pf);
  const auto struc = LoadStructure(r, s);
  constexpr auto linear_multiloop = true;  // Use same efn calculation as DP.
  efn2(data_.get(), struc.get(), 1, linear_multiloop, static_cast<std::ostream*>(nullptr));
  // TODO(2): convert ctds and structure?
  // Note: build_structure is ignored - RNAstructure bridge doesn't support Structure output.
  return {ToEnergy(struc->GetEnergy(1)), Ctds(), nullptr};
}

FoldResult RNAstructure::Fold(const Primary& r, std::optional<MfeAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const {
  VerifySupported(alg, cfg, pf);
  VerifySupported(trace_cfg);
  dp_state_t state;
  return FoldAndDpTable(r, &state);
}

FoldResult RNAstructure::FoldAndDpTable(const Primary& r, dp_state_t* dp_state) const {
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
  return {.mfe = {.dp{}, .energy = ToEnergy(structure->GetEnergy(1))},
      .tb = trace::TraceResult(StructureToSecondary(*structure), Ctds())};
}

int64_t RNAstructure::Subopt(const Primary& r, std::optional<MfeAlg> mfe_alg,
    std::optional<SuboptAlg> alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
    const subopt::SuboptCallback& fn, subopt::SuboptCfg subopt_cfg) const {
  auto res = SuboptIntoVector(r, mfe_alg, alg, cfg, pf, subopt_cfg);
  for (const auto& subopt : res) fn(subopt);
  return static_cast<int64_t>(res.size());
}

std::vector<subopt::SuboptResult> RNAstructure::SuboptIntoVector(const Primary& r,
    std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg) const {
  VerifySupported(mfe_alg, cfg, pf);
  verify(!alg.has_value(), "Bridge for RNAstructure does not support subopt algorithm selection");
  VerifySupported(subopt_cfg);
  const auto structure = LoadStructure(r);
  // Arguments: structure, data tables, percentage delta, absolute delta, nullptr, nullptr, false
  alltrace(
      structure.get(), data_.get(), 100, FromEnergy(subopt_cfg.delta), nullptr, nullptr, false);
  return StructureToSuboptVector(*structure);
}

pfn::PfnResult RNAstructure::Pfn(const Primary& r, std::optional<PfnAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf) const {
  VerifySupported(alg, cfg, pf);
  const auto structure = LoadStructure(r);
  auto state = RunPfn(structure.get(), data_.get());
  const int N = r.size();

  // RNAstructure partition values are stored in natural log space.
  auto p = BoltzSums(N, 0);
  auto q = BoltzEnergy(exp(state.w5[N]));
  for (int i = 1; i <= N; ++i) {
    for (int j = i; j < N + i; ++j) {
      const int adjusted = j > N ? j - N - 1 : j - 1;
      p[i - 1][adjusted] = BoltzEnergy(exp(state.v.f(i, j)));
    }
  }

  auto prob = BoltzProbs(N, 0);
  for (int i = 0; i < N; ++i) {
    for (int j = i; j < N; ++j) {
      prob[i][j] = BoltzEnergy(calculateprobability(i + 1, j + 1, &state.v, state.w5.get(),
          structure.get(), state.pfdata.get(), state.lfce.get(), state.mod.get(), state.scaling,
          state.fce.get()));
    }
  }
  return {.state{}, .pfn{std::move(p), q, std::move(prob)}};
}

std::vector<subopt::SuboptResult> RNAstructure::StochasticSampleIntoVector(
    const Primary& r, int num_samples) const {
  const auto structure = LoadStructure(r);
  auto state = RunPfn(structure.get(), data_.get());
  stochastictraceback(&state.w, &state.wmb, &state.wmbl, &state.wcoax, &state.wl, &state.wlc,
      &state.v, state.fce.get(), state.w3.get(), state.w5.get(), state.scaling, state.lfce.get(),
      state.mod.get(), state.pfdata.get(), num_samples, structure.get());
  return StructureToSuboptVector(*structure);
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const Primary& r) const {
  auto struc = std::make_unique<structure>();
  struc->SetThermodynamicDataTable(data_.get());
  struc->SetSequence(r.ToSeq());
  verify(struc->GetSequenceLength() == static_cast<int>(r.size()), "BUG: structure not loaded");
  return struc;
}

std::unique_ptr<structure> RNAstructure::LoadStructure(const Primary& r, const Secondary& s) const {
  auto struc = LoadStructure(r);
  struc->AddStructure();
  for (int i = 0; i < s.size(); ++i) {
    if (!s.IsOpeningPair(i)) continue;
    struc->SetPair(i + 1, s[i] + 1);
  }
  return struc;
}

RNAstructure RNAstructure::FromArgParse(const ArgParse& args) {
  return {args.Get(OPT_RNASTRUCTURE_DATA), /*use_lyngso=*/false};
}

}  // namespace mrna::bridge
