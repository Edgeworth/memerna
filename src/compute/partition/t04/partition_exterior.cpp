// Copyright 2016 Eliot Courtney.
#include <memory>

#include "compute/boltz_dp.h"
#include "compute/energy/t04/model.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::part::t04 {

BoltzExtArray PartitionExterior(
    const Primary& r, const erg::t04::ModelPtr& em, const BoltzDpArray& dp) {
  const int N = static_cast<int>(r.size());
  auto ext = BoltzExtArray(r.size() + 1, 0);

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
      const BoltzEnergy base00 = dp[st][en][PT_P] * em->AuGuPenalty(stb, enb).Boltz();
      const BoltzEnergy base01 = dp[st][en - 1][PT_P] * em->AuGuPenalty(stb, en1b).Boltz();
      const BoltzEnergy base10 = dp[st + 1][en][PT_P] * em->AuGuPenalty(st1b, enb).Boltz();
      const BoltzEnergy base11 = dp[st + 1][en - 1][PT_P] * em->AuGuPenalty(st1b, en1b).Boltz();

      // (   )<   >
      BoltzEnergy val = base00 * ext[en + 1][PTEXT_R];
      ext[st][PTEXT_R] += val;
      if (IsGuPair(stb, enb))
        ext[st][PTEXT_R_GU] += val;
      else
        ext[st][PTEXT_R_WC] += val;

      // (   )3<   > 3'
      ext[st][PTEXT_R] += base01 * em->dangle3[en1b][enb][stb].Boltz() * ext[en + 1][PTEXT_R];
      // 5(   )<   > 5'
      ext[st][PTEXT_R] += base10 * em->dangle5[enb][stb][st1b].Boltz() * ext[en + 1][PTEXT_R];
      // .(   ).<   > Terminal mismatch
      ext[st][PTEXT_R] +=
          base11 * em->terminal[en1b][enb][stb][st1b].Boltz() * ext[en + 1][PTEXT_R];
      // .(   ).<(   ) > Left coax
      val = base11 * em->MismatchCoaxial(en1b, enb, stb, st1b).Boltz();
      ext[st][PTEXT_R] += val * ext[en + 1][PTEXT_R_GU];
      ext[st][PTEXT_R] += val * ext[en + 1][PTEXT_R_WC];

      // (   )<.(   ). > Right coax forward
      ext[st][PTEXT_R] += base00 * ext[en + 1][PT_U_RC];
      // (   )<.( * ). > Right coax backward
      ext[st][PT_U_RC] +=
          base11 * em->MismatchCoaxial(en1b, enb, stb, st1b).Boltz() * ext[en + 1][PTEXT_R];

      // (   )(<   ) > Flush coax
      ext[st][PTEXT_R] +=
          base01 * em->stack[en1b][enb][WcPair(enb)][stb].Boltz() * ext[en][PTEXT_R_WC];
      if (IsGu(enb))
        ext[st][PTEXT_R] +=
            base01 * em->stack[en1b][enb][GuPair(enb)][stb].Boltz() * ext[en][PTEXT_R_GU];
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
      const BoltzEnergy base00 = dp[st][en][PT_P] * em->AuGuPenalty(stb, enb).Boltz();
      const BoltzEnergy base01 = dp[st][en - 1][PT_P] * em->AuGuPenalty(stb, en1b).Boltz();
      const BoltzEnergy base10 = dp[st + 1][en][PT_P] * em->AuGuPenalty(st1b, enb).Boltz();
      const BoltzEnergy base11 = dp[st + 1][en - 1][PT_P] * em->AuGuPenalty(st1b, en1b).Boltz();
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

      // <   >(   )3 3'
      ext[en][PTEXT_L] += base01 * em->dangle3[en1b][enb][stb].Boltz() * ptextl;
      // <   >5(   ) 5'
      ext[en][PTEXT_L] += base10 * em->dangle5[enb][stb][st1b].Boltz() * ptextl;
      // <   >.(   ). Terminal mismatch
      ext[en][PTEXT_L] += base11 * em->terminal[en1b][enb][stb][st1b].Boltz() * ptextl;
      // <  (   )>.(   ). Right coax
      val = base11 * em->MismatchCoaxial(en1b, enb, stb, st1b).Boltz();
      ext[en][PTEXT_L] += val * ptextlgu;
      ext[en][PTEXT_L] += val * ptextlwc;

      // <  .(   ).>(   ) Left coax forward
      ext[en][PTEXT_L] += base00 * ptextllcoaxx;
      // <  .( * ).>(   ) Left coax backward
      ext[en][PTEXT_L_LCOAX] += base11 * em->MismatchCoaxial(en1b, enb, stb, st1b).Boltz() * ptextl;

      // < (   >)(   ) Flush coax
      ext[en][PTEXT_L] +=
          base10 * em->stack[stb][st1b][enb][WcPair(stb)].Boltz() * ext[st][PTEXT_L_WC];
      if (IsGu(stb))
        ext[en][PTEXT_L] +=
            base10 * em->stack[stb][st1b][enb][GuPair(stb)].Boltz() * ext[st][PTEXT_L_GU];
    }
  }

  return ext;
}

}  // namespace mrna::part::t04
