// Copyright 2026 Eliot Courtney.
#include "backends/canon/pfn/pfn_exterior.h"

#include "api/energy/energy_cfg.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
void PfnExterior(const Primary& r, const Model& m, PfnState& state, const PseudofreeCfgT& pf) {
  const int N = r.size();
  pf.Verify(r);

  const auto& dp = state.dp;
  state.ext = BoltzExtArray(r.size() + 1, ZERO_B);
  auto& ext = state.ext;

  ext[N][PTEXT_R] = ONE_B;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    ext[st][PTEXT_R] += ext[st + 1][PTEXT_R] * pf.Unpaired(st).Boltz();
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const BoltzEnergy base00 = dp[st][en][PT_P] * m.AuGuPenalty(stb, enb).Boltz();
      const BoltzEnergy base01 = dp[st][en - 1][PT_P] * m.AuGuPenalty(stb, en1b).Boltz();
      const BoltzEnergy base10 = dp[st + 1][en][PT_P] * m.AuGuPenalty(st1b, enb).Boltz();
      const BoltzEnergy base11 = dp[st + 1][en - 1][PT_P] * m.AuGuPenalty(st1b, en1b).Boltz();

      // (   )<   >
      BoltzEnergy val = base00 * ext[en + 1][PTEXT_R];

      if constexpr (CFG.UseD2()) {
        if (st != 0 && en != N - 1) {
          // (   )<   > Terminal mismatch
          val *= m.terminal[enb][r[en + 1]][r[st - 1]][stb].Boltz();
        } else if (en != N - 1) {
          // (   )<3   > 3'
          val *= m.dangle3[enb][r[en + 1]][stb].Boltz();
        } else if (st != 0) {
          // 5(   )<   > 5'
          val *= m.dangle5[enb][r[st - 1]][stb].Boltz();
        }
      }

      ext[st][PTEXT_R] += val;

      if constexpr (CFG.UseDangleMismatch()) {
        // (   )3<   > 3'
        ext[st][PTEXT_R] +=
            base01 * (m.dangle3[en1b][enb][stb] + pf.Unpaired(en)).Boltz() * ext[en + 1][PTEXT_R];
        // 5(   )<   > 5'
        ext[st][PTEXT_R] +=
            base10 * (m.dangle5[enb][stb][st1b] + pf.Unpaired(st)).Boltz() * ext[en + 1][PTEXT_R];
        // .(   ).<   > Terminal mismatch
        ext[st][PTEXT_R] += base11 *
            (m.terminal[en1b][enb][stb][st1b] + pf.Unpaired(st) + pf.Unpaired(en)).Boltz() *
            ext[en + 1][PTEXT_R];
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // .(   ).<(   ) > Left coax
        for (int piv = en + HAIRPIN_MIN_SZ + 2; piv < N; ++piv) {
          const Base prb = r[en + 1];
          const Base pb = r[piv];
          ext[st][PTEXT_R] += base11 * dp[en + 1][piv][PT_P] * m.AuGuPenalty(prb, pb).Boltz() *
              (m.MismatchCoaxial(en1b, enb, stb, st1b) + pf.Unpaired(st) + pf.Unpaired(en))
                  .Boltz() *
              ext[piv + 1][PTEXT_R];
        }

        // (   )<.(   ). > Right coax forward
        for (int piv = en + HAIRPIN_MIN_SZ + 4; piv < N; ++piv) {
          const Base prb = r[en + 1];
          const Base pr1b = r[en + 2];
          const Base pl1b = r[piv - 1];
          const Base pb = r[piv];
          ext[st][PTEXT_R] += base00 * dp[en + 2][piv - 1][PT_P] *
              m.AuGuPenalty(pr1b, pl1b).Boltz() *
              (m.MismatchCoaxial(pl1b, pb, prb, pr1b) + pf.Unpaired(en + 1) + pf.Unpaired(piv))
                  .Boltz() *
              ext[piv + 1][PTEXT_R];
        }

        // (   )(<   ) > Flush coax
        for (int piv = en + HAIRPIN_MIN_SZ + 1; piv < N; ++piv) {
          const Base pb = r[piv];
          ext[st][PTEXT_R] += base01 * dp[en][piv][PT_P] * m.AuGuPenalty(enb, pb).Boltz() *
              m.stack[en1b][enb][pb][stb].Boltz() * ext[piv + 1][PTEXT_R];
        }
      }
    }
  }

  ext[0][PTEXT_L] = pf.Unpaired(0).Boltz();
  for (int en = 1; en < N; ++en) {
    // Case: No pair ending here
    ext[en][PTEXT_L] += ext[en - 1][PTEXT_L] * pf.Unpaired(en).Boltz();
    for (int st = 0; st < en - HAIRPIN_MIN_SZ; ++st) {
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const BoltzEnergy base00 = dp[st][en][PT_P] * m.AuGuPenalty(stb, enb).Boltz();
      const BoltzEnergy base01 = dp[st][en - 1][PT_P] * m.AuGuPenalty(stb, en1b).Boltz();
      const BoltzEnergy base10 = dp[st + 1][en][PT_P] * m.AuGuPenalty(st1b, enb).Boltz();
      const BoltzEnergy base11 = dp[st + 1][en - 1][PT_P] * m.AuGuPenalty(st1b, en1b).Boltz();

      const BoltzEnergy ptextl = st ? ext[st - 1][PTEXT_L] : BoltzEnergy{1};

      // <   >(   )
      BoltzEnergy val = base00 * ptextl;

      if constexpr (CFG.UseD2()) {
        if (st != 0 && en != N - 1) {
          // <   m>(   )m Terminal mismatch
          val *= m.terminal[enb][r[en + 1]][r[st - 1]][stb].Boltz();
        } else if (en != N - 1) {
          // <   >(   )3 3'
          val *= m.dangle3[enb][r[en + 1]][stb].Boltz();
        } else if (st != 0) {
          // <   5>(   ) 5'
          val *= m.dangle5[enb][r[st - 1]][stb].Boltz();
        }
      }
      ext[en][PTEXT_L] += val;

      if constexpr (CFG.UseDangleMismatch()) {
        // <   >(   )3 3'
        ext[en][PTEXT_L] += base01 * (m.dangle3[en1b][enb][stb] + pf.Unpaired(en)).Boltz() * ptextl;
        // <   >5(   ) 5'
        ext[en][PTEXT_L] += base10 * (m.dangle5[enb][stb][st1b] + pf.Unpaired(st)).Boltz() * ptextl;
        // <   >.(   ). Terminal mismatch
        ext[en][PTEXT_L] += base11 *
            (m.terminal[en1b][enb][stb][st1b] + pf.Unpaired(st) + pf.Unpaired(en)).Boltz() * ptextl;
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // <  (   )>.(   ). Right coax
        for (int piv = 0; piv <= st - HAIRPIN_MIN_SZ - 2; ++piv) {
          const BoltzEnergy prefix = piv ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
          const Base pb = r[piv];
          const Base pl1b = r[st - 1];
          ext[en][PTEXT_L] += prefix * dp[piv][st - 1][PT_P] * m.AuGuPenalty(pb, pl1b).Boltz() *
              base11 *
              (m.MismatchCoaxial(en1b, enb, stb, st1b) + pf.Unpaired(st) + pf.Unpaired(en)).Boltz();
        }

        // <  .(   ).>(   ) Left coax forward
        for (int piv = 0; piv <= st - HAIRPIN_MIN_SZ - 4; ++piv) {
          const BoltzEnergy prefix = piv ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
          const Base pb = r[piv];
          const Base prb = r[piv + 1];
          const Base pl1b = r[st - 2];
          const Base plb = r[st - 1];
          ext[en][PTEXT_L] += prefix * dp[piv + 1][st - 2][PT_P] *
              m.AuGuPenalty(prb, pl1b).Boltz() * base00 *
              (m.MismatchCoaxial(pl1b, plb, pb, prb) + pf.Unpaired(piv) + pf.Unpaired(st - 1))
                  .Boltz();
        }

        // < (   >)(   ) Flush coax
        for (int piv = 0; piv <= st - HAIRPIN_MIN_SZ - 1; ++piv) {
          const BoltzEnergy prefix = piv ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
          const Base pb = r[piv];
          ext[en][PTEXT_L] += prefix * dp[piv][st][PT_P] * m.AuGuPenalty(pb, stb).Boltz() * base10 *
              m.stack[stb][st1b][enb][pb].Boltz();
        }
      }
    }
  }
}

#define MRNA_INSTANTIATE(CFG, PF) \
  template void PfnExterior<CFG, PF>(const Primary&, const Model&, PfnState&, const PF&);
MRNA_STATIC_FOR_EACH_PFN_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::canon
