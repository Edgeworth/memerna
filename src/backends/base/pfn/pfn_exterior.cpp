// Copyright 2016 E.
#include "api/energy/energy_cfg.h"
#include "backends/base/pfn/pfn.h"
#include "backends/common/base/model_base.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna::md::base {

void PfnExterior(const Primary& r, const ModelBase& m, PfnState& state) {
  const int N = static_cast<int>(r.size());

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false},  // Bulge states with partition function doesn't make sense.
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m.cfg());

  const auto& dp = state.dp;
  state.ext = BoltzExtArray(r.size() + 1, 0);
  auto& ext = state.ext;

  ext[N][PTEXT_R] = 1;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    ext[st][PTEXT_R] += ext[st + 1][PTEXT_R];
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
      ext[st][PTEXT_R] += val;
      if (IsGuPair(stb, enb))
        ext[st][PTEXT_R_GU] += val;
      else
        ext[st][PTEXT_R_WC] += val;

      if (m.cfg().UseDangleMismatch()) {
        // (   )3<   > 3'
        ext[st][PTEXT_R] += base01 * m.dangle3[en1b][enb][stb].Boltz() * ext[en + 1][PTEXT_R];
        // 5(   )<   > 5'
        ext[st][PTEXT_R] += base10 * m.dangle5[enb][stb][st1b].Boltz() * ext[en + 1][PTEXT_R];
        // .(   ).<   > Terminal mismatch
        ext[st][PTEXT_R] +=
            base11 * m.terminal[en1b][enb][stb][st1b].Boltz() * ext[en + 1][PTEXT_R];
      }

      if (m.cfg().UseCoaxialStacking()) {
        // .(   ).<(   ) > Left coax
        val = base11 * m.MismatchCoaxial(en1b, enb, stb, st1b).Boltz();
        ext[st][PTEXT_R] += val * ext[en + 1][PTEXT_R_GU];
        ext[st][PTEXT_R] += val * ext[en + 1][PTEXT_R_WC];

        // (   )<.(   ). > Right coax forward
        ext[st][PTEXT_R] += base00 * ext[en + 1][PTEXT_R_RC];
        // (   )<.( * ). > Right coax backward
        ext[st][PTEXT_R_RC] +=
            base11 * m.MismatchCoaxial(en1b, enb, stb, st1b).Boltz() * ext[en + 1][PTEXT_R];

        // (   )(<   ) > Flush coax
        ext[st][PTEXT_R] +=
            base01 * m.stack[en1b][enb][WcPair(enb)][stb].Boltz() * ext[en][PTEXT_R_WC];
        if (IsGu(enb))
          ext[st][PTEXT_R] +=
              base01 * m.stack[en1b][enb][GuPair(enb)][stb].Boltz() * ext[en][PTEXT_R_GU];
      }
    }
  }

  ext[0][PTEXT_L] = 1;
  for (int en = 1; en < N; ++en) {
    // Case: No pair ending here
    ext[en][PTEXT_L] += ext[en - 1][PTEXT_L];
    for (int st = 0; st < en - HAIRPIN_MIN_SZ; ++st) {
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const BoltzEnergy base00 = dp[st][en][PT_P] * m.AuGuPenalty(stb, enb).Boltz();
      const BoltzEnergy base01 = dp[st][en - 1][PT_P] * m.AuGuPenalty(stb, en1b).Boltz();
      const BoltzEnergy base10 = dp[st + 1][en][PT_P] * m.AuGuPenalty(st1b, enb).Boltz();
      const BoltzEnergy base11 = dp[st + 1][en - 1][PT_P] * m.AuGuPenalty(st1b, en1b).Boltz();
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
      ext[en][PTEXT_L] += val;
      if (IsGuPair(stb, enb))
        ext[en][PTEXT_L_GU] += val;
      else
        ext[en][PTEXT_L_WC] += val;

      if (m.cfg().UseDangleMismatch()) {
        // <   >(   )3 3'
        ext[en][PTEXT_L] += base01 * m.dangle3[en1b][enb][stb].Boltz() * ptextl;
        // <   >5(   ) 5'
        ext[en][PTEXT_L] += base10 * m.dangle5[enb][stb][st1b].Boltz() * ptextl;
        // <   >.(   ). Terminal mismatch
        ext[en][PTEXT_L] += base11 * m.terminal[en1b][enb][stb][st1b].Boltz() * ptextl;
      }

      if (m.cfg().UseCoaxialStacking()) {
        // <  (   )>.(   ). Right coax
        val = base11 * m.MismatchCoaxial(en1b, enb, stb, st1b).Boltz();
        ext[en][PTEXT_L] += val * ptextlgu;
        ext[en][PTEXT_L] += val * ptextlwc;

        // <  .(   ).>(   ) Left coax forward
        ext[en][PTEXT_L] += base00 * ptextllcoaxx;
        // <  .( * ).>(   ) Left coax backward
        ext[en][PTEXT_L_LCOAX] += base11 * m.MismatchCoaxial(en1b, enb, stb, st1b).Boltz() * ptextl;

        // < (   >)(   ) Flush coax
        ext[en][PTEXT_L] +=
            base10 * m.stack[stb][st1b][enb][WcPair(stb)].Boltz() * ext[st][PTEXT_L_WC];
        if (IsGu(stb))
          ext[en][PTEXT_L] +=
              base10 * m.stack[stb][st1b][enb][GuPair(stb)].Boltz() * ext[st][PTEXT_L_GU];
      }
    }
  }
}

}  // namespace mrna::md::base
