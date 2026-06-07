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
void PfnExterior(
    const Primary& r, const BoltzModel& bm, PfnState& state, const PseudofreeCfgT& pf) {
  const int N = r.size();
  pf.Verify(r);
  const auto bpf = Boltz(pf);

  const auto& dp = state.dp;
  state.ext = BoltzExtArray(r.size() + 1, ZERO_B);
  auto& ext = state.ext;

  ext[N][PTEXT_R] = ONE_B;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    ext[st][PTEXT_R] += ext[st + 1][PTEXT_R] * bpf.Unpaired(st);
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const BoltzEnergy base00 = dp[st][en][PT_P] * bm.AuGuPenalty(stb, enb);
      const BoltzEnergy base01 = dp[st][en - 1][PT_P] * bm.AuGuPenalty(stb, en1b);
      const BoltzEnergy base10 = dp[st + 1][en][PT_P] * bm.AuGuPenalty(st1b, enb);
      const BoltzEnergy base11 = dp[st + 1][en - 1][PT_P] * bm.AuGuPenalty(st1b, en1b);

      // (   )<   >
      BoltzEnergy val = base00 * ext[en + 1][PTEXT_R];

      if constexpr (CFG.UseD2()) {
        if (st != 0 && en != N - 1) {
          // (   )<   > Terminal mismatch
          val *= bm.terminal[enb][r[en + 1]][r[st - 1]][stb];
        } else if (en != N - 1) {
          // (   )<3   > 3'
          val *= bm.dangle3[enb][r[en + 1]][stb];
        } else if (st != 0) {
          // 5(   )<   > 5'
          val *= bm.dangle5[enb][r[st - 1]][stb];
        }
      }

      ext[st][PTEXT_R] += val;

      if constexpr (CFG.UseDangleMismatch()) {
        // (   )3<   > 3'
        ext[st][PTEXT_R] +=
            base01 * bm.dangle3[en1b][enb][stb] * bpf.Unpaired(en) * ext[en + 1][PTEXT_R];
        // 5(   )<   > 5'
        ext[st][PTEXT_R] +=
            base10 * bm.dangle5[enb][stb][st1b] * bpf.Unpaired(st) * ext[en + 1][PTEXT_R];
        // .(   ).<   > Terminal mismatch
        ext[st][PTEXT_R] += base11 * bm.terminal[en1b][enb][stb][st1b] *
            bpf.Unpaired(st) * bpf.Unpaired(en) * ext[en + 1][PTEXT_R];
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // .(   ).<(   ) > Left coax
        for (int piv = en + HAIRPIN_MIN_SZ + 2; piv < N; ++piv) {
          const Base prb = r[en + 1];
          const Base pb = r[piv];
          ext[st][PTEXT_R] += base11 * dp[en + 1][piv][PT_P] * bm.AuGuPenalty(prb, pb) *
              bm.MismatchCoaxial(en1b, enb, stb, st1b) * bpf.Unpaired(st) *
              bpf.Unpaired(en) *
              ext[piv + 1][PTEXT_R];
        }

        // (   )<.(   ). > Right coax forward
        for (int piv = en + HAIRPIN_MIN_SZ + 4; piv < N; ++piv) {
          const Base prb = r[en + 1];
          const Base pr1b = r[en + 2];
          const Base pl1b = r[piv - 1];
          const Base pb = r[piv];
          ext[st][PTEXT_R] += base00 * dp[en + 2][piv - 1][PT_P] *
              bm.AuGuPenalty(pr1b, pl1b) * bm.MismatchCoaxial(pl1b, pb, prb, pr1b) *
              bpf.Unpaired(en + 1) * bpf.Unpaired(piv) *
              ext[piv + 1][PTEXT_R];
        }

        // (   )(<   ) > Flush coax
        for (int piv = en + HAIRPIN_MIN_SZ + 1; piv < N; ++piv) {
          const Base pb = r[piv];
          ext[st][PTEXT_R] += base01 * dp[en][piv][PT_P] * bm.AuGuPenalty(enb, pb) *
              bm.stack[en1b][enb][pb][stb] * ext[piv + 1][PTEXT_R];
        }
      }
    }
  }

  ext[0][PTEXT_L] = bpf.Unpaired(0);
  for (int en = 1; en < N; ++en) {
    // Case: No pair ending here
    ext[en][PTEXT_L] += ext[en - 1][PTEXT_L] * bpf.Unpaired(en);
    for (int st = 0; st < en - HAIRPIN_MIN_SZ; ++st) {
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const BoltzEnergy base00 = dp[st][en][PT_P] * bm.AuGuPenalty(stb, enb);
      const BoltzEnergy base01 = dp[st][en - 1][PT_P] * bm.AuGuPenalty(stb, en1b);
      const BoltzEnergy base10 = dp[st + 1][en][PT_P] * bm.AuGuPenalty(st1b, enb);
      const BoltzEnergy base11 = dp[st + 1][en - 1][PT_P] * bm.AuGuPenalty(st1b, en1b);

      const BoltzEnergy ptextl = st ? ext[st - 1][PTEXT_L] : BoltzEnergy{1};

      // <   >(   )
      BoltzEnergy val = base00 * ptextl;

      if constexpr (CFG.UseD2()) {
        if (st != 0 && en != N - 1) {
          // <   m>(   )m Terminal mismatch
          val *= bm.terminal[enb][r[en + 1]][r[st - 1]][stb];
        } else if (en != N - 1) {
          // <   >(   )3 3'
          val *= bm.dangle3[enb][r[en + 1]][stb];
        } else if (st != 0) {
          // <   5>(   ) 5'
          val *= bm.dangle5[enb][r[st - 1]][stb];
        }
      }
      ext[en][PTEXT_L] += val;

      if constexpr (CFG.UseDangleMismatch()) {
        // <   >(   )3 3'
        ext[en][PTEXT_L] += base01 * bm.dangle3[en1b][enb][stb] * bpf.Unpaired(en) * ptextl;
        // <   >5(   ) 5'
        ext[en][PTEXT_L] += base10 * bm.dangle5[enb][stb][st1b] * bpf.Unpaired(st) * ptextl;
        // <   >.(   ). Terminal mismatch
        ext[en][PTEXT_L] += base11 * bm.terminal[en1b][enb][stb][st1b] *
            bpf.Unpaired(st) * bpf.Unpaired(en) * ptextl;
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // <  (   )>.(   ). Right coax
        for (int piv = 0; piv <= st - HAIRPIN_MIN_SZ - 2; ++piv) {
          const BoltzEnergy prefix = piv ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
          const Base pb = r[piv];
          const Base pl1b = r[st - 1];
          ext[en][PTEXT_L] += prefix * dp[piv][st - 1][PT_P] * bm.AuGuPenalty(pb, pl1b) *
              base11 * bm.MismatchCoaxial(en1b, enb, stb, st1b) * bpf.Unpaired(st) *
              bpf.Unpaired(en);
        }

        // <  .(   ).>(   ) Left coax forward
        for (int piv = 0; piv <= st - HAIRPIN_MIN_SZ - 4; ++piv) {
          const BoltzEnergy prefix = piv ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
          const Base pb = r[piv];
          const Base prb = r[piv + 1];
          const Base pl1b = r[st - 2];
          const Base plb = r[st - 1];
          ext[en][PTEXT_L] += prefix * dp[piv + 1][st - 2][PT_P] *
              bm.AuGuPenalty(prb, pl1b) * base00 * bm.MismatchCoaxial(pl1b, plb, pb, prb) *
              bpf.Unpaired(piv) * bpf.Unpaired(st - 1);
        }

        // < (   >)(   ) Flush coax
        for (int piv = 0; piv <= st - HAIRPIN_MIN_SZ - 1; ++piv) {
          const BoltzEnergy prefix = piv ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
          const Base pb = r[piv];
          ext[en][PTEXT_L] += prefix * dp[piv][st][PT_P] * bm.AuGuPenalty(pb, stb) * base10 *
              bm.stack[stb][st1b][enb][pb];
        }
      }
    }
  }
}

#define MRNA_INSTANTIATE(CFG, PF) \
  template void PfnExterior<CFG, PF>(const Primary&, const BoltzModel&, PfnState&, const PF&);
MRNA_STATIC_FOR_EACH_PFN_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::canon
