// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_COMMON_STATIC_H_
#define BACKENDS_COMMON_STATIC_H_

#include <utility>

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "util/error.h"

#define MRNA_STATIC_CFG(LONELY_PAIRS, BULGE_STATES, CTD)                                     \
  (::mrna::erg::EnergyCfg{.lonely_pairs = ::mrna::erg::EnergyCfg::LonelyPairs::LONELY_PAIRS, \
      .bulge_states = BULGE_STATES,                                                          \
      .ctd = ::mrna::erg::EnergyCfg::Ctd::CTD})

#define MRNA_STATIC_FOR_EACH_PF(M, CFG) \
  M(CFG, ::mrna::erg::NoPseudofreeCfg)  \
  M(CFG, ::mrna::erg::PseudofreeCfg)

#define MRNA_STATIC_FOR_EACH_CTD_PF(M, LONELY_PAIRS, BULGE_STATES)                 \
  MRNA_STATIC_FOR_EACH_PF(M, MRNA_STATIC_CFG(LONELY_PAIRS, BULGE_STATES, NONE))    \
  MRNA_STATIC_FOR_EACH_PF(M, MRNA_STATIC_CFG(LONELY_PAIRS, BULGE_STATES, D2))      \
  MRNA_STATIC_FOR_EACH_PF(M, MRNA_STATIC_CFG(LONELY_PAIRS, BULGE_STATES, NO_COAX)) \
  MRNA_STATIC_FOR_EACH_PF(M, MRNA_STATIC_CFG(LONELY_PAIRS, BULGE_STATES, ALL))

#define MRNA_STATIC_FOR_EACH_BULGE_PF(M, LONELY_PAIRS) \
  MRNA_STATIC_FOR_EACH_CTD_PF(M, LONELY_PAIRS, false)  \
  MRNA_STATIC_FOR_EACH_CTD_PF(M, LONELY_PAIRS, true)

#define MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(M) \
  MRNA_STATIC_FOR_EACH_BULGE_PF(M, OFF)       \
  MRNA_STATIC_FOR_EACH_BULGE_PF(M, HEURISTIC) \
  MRNA_STATIC_FOR_EACH_BULGE_PF(M, ON)

#define MRNA_STATIC_FOR_EACH_PFN_ENERGY_CFG_PF(M)  \
  MRNA_STATIC_FOR_EACH_CTD_PF(M, OFF, false)       \
  MRNA_STATIC_FOR_EACH_CTD_PF(M, HEURISTIC, false) \
  MRNA_STATIC_FOR_EACH_CTD_PF(M, ON, false)

namespace mrna::md {

template <erg::EnergyCfg::LonelyPairs LonelyPairs, bool BulgeStates, erg::EnergyCfg::Ctd Ctd,
    class Fn>
decltype(auto) Dispatch(Fn&& fn) {
  constexpr erg::EnergyCfg CFG{
      .lonely_pairs = LonelyPairs,
      .bulge_states = BulgeStates,
      .ctd = Ctd,
  };
  return std::forward<Fn>(fn).template operator()<CFG>();
}

template <erg::EnergyCfg::LonelyPairs LonelyPairs, bool BulgeStates, class Fn>
decltype(auto) Dispatch(erg::EnergyCfg::Ctd ctd, Fn&& fn) {
  switch (ctd) {
  case erg::EnergyCfg::Ctd::NONE:
    return Dispatch<LonelyPairs, BulgeStates, erg::EnergyCfg::Ctd::NONE>(std::forward<Fn>(fn));
  case erg::EnergyCfg::Ctd::D2:
    return Dispatch<LonelyPairs, BulgeStates, erg::EnergyCfg::Ctd::D2>(std::forward<Fn>(fn));
  case erg::EnergyCfg::Ctd::NO_COAX:
    return Dispatch<LonelyPairs, BulgeStates, erg::EnergyCfg::Ctd::NO_COAX>(std::forward<Fn>(fn));
  case erg::EnergyCfg::Ctd::ALL:
    return Dispatch<LonelyPairs, BulgeStates, erg::EnergyCfg::Ctd::ALL>(std::forward<Fn>(fn));
  }
  unreachable();
}

template <erg::EnergyCfg::LonelyPairs LonelyPairs, class Fn>
decltype(auto) Dispatch(bool bulge_states, erg::EnergyCfg::Ctd ctd, Fn&& fn) {
  if (bulge_states) return Dispatch<LonelyPairs, true>(ctd, std::forward<Fn>(fn));
  return Dispatch<LonelyPairs, false>(ctd, std::forward<Fn>(fn));
}

template <class Fn>
decltype(auto) Dispatch(erg::EnergyCfg cfg, Fn&& fn) {
  switch (cfg.lonely_pairs) {
  case erg::EnergyCfg::LonelyPairs::OFF:
    return Dispatch<erg::EnergyCfg::LonelyPairs::OFF>(
        cfg.bulge_states, cfg.ctd, std::forward<Fn>(fn));
  case erg::EnergyCfg::LonelyPairs::HEURISTIC:
    return Dispatch<erg::EnergyCfg::LonelyPairs::HEURISTIC>(
        cfg.bulge_states, cfg.ctd, std::forward<Fn>(fn));
  case erg::EnergyCfg::LonelyPairs::ON:
    return Dispatch<erg::EnergyCfg::LonelyPairs::ON>(
        cfg.bulge_states, cfg.ctd, std::forward<Fn>(fn));
  }
  unreachable();
}

template <class Fn>
decltype(auto) Dispatch(erg::EnergyCfgPfn cfg, Fn&& fn) {
  switch (cfg.lonely_pairs) {
  case erg::EnergyCfg::LonelyPairs::OFF:
    return Dispatch<erg::EnergyCfg::LonelyPairs::OFF, false>(cfg.ctd, std::forward<Fn>(fn));
  case erg::EnergyCfg::LonelyPairs::HEURISTIC:
    return Dispatch<erg::EnergyCfg::LonelyPairs::HEURISTIC, false>(cfg.ctd, std::forward<Fn>(fn));
  case erg::EnergyCfg::LonelyPairs::ON:
    return Dispatch<erg::EnergyCfg::LonelyPairs::ON, false>(cfg.ctd, std::forward<Fn>(fn));
  }
  unreachable();
}

template <class Fn>
decltype(auto) Dispatch(const erg::PseudofreeCfg& pf, Fn&& fn) {
  if (pf.Empty())
    return std::forward<Fn>(fn).template operator()<erg::NoPseudofreeCfg>(erg::NoPseudofreeCfg{pf});
  return std::forward<Fn>(fn).template operator()<erg::PseudofreeCfg>(pf);
}

template <class EnergyCfgT, class Fn>
decltype(auto) Dispatch(EnergyCfgT cfg, const erg::PseudofreeCfg& pf, Fn&& fn) {
  return Dispatch(cfg, [&]<erg::EnergyCfg CFG>() -> decltype(auto) {
    return Dispatch(
        pf, [&]<class PseudofreeCfgT>(const PseudofreeCfgT& static_pf) -> decltype(auto) {
          return std::forward<Fn>(fn).template operator()<CFG, PseudofreeCfgT>(static_pf);
        });
  });
}

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_STATIC_H_
