// Copyright 2016 Eliot Courtney.
#include "backends/base/pfn/pfn_exterior.h"

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/dp.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

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
          // (   )<   > Terminal mismatch - U
          val *= bm.terminal[enb][r[en + 1]][r[st - 1]][stb];
        } else if (en != N - 1) {
          // (   )<3   > 3' - U
          val *= bm.dangle3[enb][r[en + 1]][stb];
        } else if (st != 0) {
          // 5(   )<   > 5' - U
          val *= bm.dangle5[enb][r[st - 1]][stb];
        }
      }

      ext[st][PTEXT_R] += val;
      if (IsGuPair(stb, enb))
        ext[st][PTEXT_R_GU] += val;
      else
        ext[st][PTEXT_R_WC] += val;

      if constexpr (CFG.UseDangleMismatch()) {
        // (   )3<   > 3'
        ext[st][PTEXT_R] +=
            base01 * bm.dangle3[en1b][enb][stb] * bpf.Unpaired(en) * ext[en + 1][PTEXT_R];
        // 5(   )<   > 5'
        ext[st][PTEXT_R] +=
            base10 * bm.dangle5[enb][stb][st1b] * bpf.Unpaired(st) * ext[en + 1][PTEXT_R];
        // .(   ).<   > Terminal mismatch
        ext[st][PTEXT_R] += base11 * bm.terminal[en1b][enb][stb][st1b] * bpf.Unpaired(st) *
            bpf.Unpaired(en) * ext[en + 1][PTEXT_R];
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // .(   ).<(   ) > Left coax
        val =
            base11 * bm.MismatchCoaxial(en1b, enb, stb, st1b) * bpf.Unpaired(st) * bpf.Unpaired(en);
        ext[st][PTEXT_R] += val * ext[en + 1][PTEXT_R_GU];
        ext[st][PTEXT_R] += val * ext[en + 1][PTEXT_R_WC];

        // (   )<.(   ). > Right coax forward
        ext[st][PTEXT_R] += base00 * ext[en + 1][PTEXT_R_RC];
        // (   )<.( * ). > Right coax backward
        ext[st][PTEXT_R_RC] += base11 * bm.MismatchCoaxial(en1b, enb, stb, st1b) *
            bpf.Unpaired(st) * bpf.Unpaired(en) * ext[en + 1][PTEXT_R];

        // (   )(<   ) > Flush coax
        ext[st][PTEXT_R] += base01 * bm.stack[en1b][enb][WcPair(enb)][stb] * ext[en][PTEXT_R_WC];
        if (IsGu(enb))
          ext[st][PTEXT_R] += base01 * bm.stack[en1b][enb][GuPair(enb)][stb] * ext[en][PTEXT_R_GU];
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
      BoltzEnergy ptextl{1};
      BoltzEnergy ptextlgu{0};
      BoltzEnergy ptextlwc{0};
      BoltzEnergy ptextllcoaxx{0};

      if (st) {
        ptextl = ext[st - 1][PTEXT_L];
        ptextlgu = ext[st - 1][PTEXT_L_GU];
        ptextlwc = ext[st - 1][PTEXT_L_WC];
        ptextllcoaxx = ext[st - 1][PTEXT_L_LCOAX];
      }

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

      if (IsGuPair(stb, enb))
        ext[en][PTEXT_L_GU] += val;
      else
        ext[en][PTEXT_L_WC] += val;

      if constexpr (CFG.UseDangleMismatch()) {
        // <   >(   )3 3'
        ext[en][PTEXT_L] += base01 * bm.dangle3[en1b][enb][stb] * bpf.Unpaired(en) * ptextl;
        // <   >5(   ) 5'
        ext[en][PTEXT_L] += base10 * bm.dangle5[enb][stb][st1b] * bpf.Unpaired(st) * ptextl;
        // <   >.(   ). Terminal mismatch
        ext[en][PTEXT_L] += base11 * bm.terminal[en1b][enb][stb][st1b] * bpf.Unpaired(st) *
            bpf.Unpaired(en) * ptextl;
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // <  (   )>.(   ). Right coax
        val =
            base11 * bm.MismatchCoaxial(en1b, enb, stb, st1b) * bpf.Unpaired(st) * bpf.Unpaired(en);
        ext[en][PTEXT_L] += val * ptextlgu;
        ext[en][PTEXT_L] += val * ptextlwc;

        // <  .(   ).>(   ) Left coax forward
        ext[en][PTEXT_L] += base00 * ptextllcoaxx;
        // <  .( * ).>(   ) Left coax backward
        ext[en][PTEXT_L_LCOAX] += base11 * bm.MismatchCoaxial(en1b, enb, stb, st1b) *
            bpf.Unpaired(st) * bpf.Unpaired(en) * ptextl;

        // < (   >)(   ) Flush coax
        ext[en][PTEXT_L] += base10 * bm.stack[stb][st1b][enb][WcPair(stb)] * ext[st][PTEXT_L_WC];
        if (IsGu(stb))
          ext[en][PTEXT_L] += base10 * bm.stack[stb][st1b][enb][GuPair(stb)] * ext[st][PTEXT_L_GU];
      }
    }
  }
}

#define MRNA_INSTANTIATE(CFG, PF) \
  template void PfnExterior<CFG, PF>(const Primary&, const BoltzModel&, PfnState&, const PF&);
MRNA_STATIC_FOR_EACH_PFN_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::base
