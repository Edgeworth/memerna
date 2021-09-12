// Copyright 2016 E.
#include "energy/energy_globals.h"
#include "parsing.h"
#include "partition/partition.h"
#include "partition/partition_globals.h"
namespace mrna {
namespace partition {
namespace internal {

using energy::Boltzmann;
using energy::gem;

void Exterior() {
  const int N = static_cast<int>(gr.size());
  gptext[N][PTEXT_R] = 1;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    gptext[st][PTEXT_R] += gptext[st + 1][PTEXT_R];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const penergy_t base00 = gpt[st][en][PT_P] * Boltzmann(gem.AuGuPenalty(stb, enb));
      const penergy_t base01 = gpt[st][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(stb, en1b));
      const penergy_t base10 = gpt[st + 1][en][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, enb));
      const penergy_t base11 = gpt[st + 1][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, en1b));

      // (   )<   >
      penergy_t val = base00 * gptext[en + 1][PTEXT_R];
      gptext[st][PTEXT_R] += val;
      if (IsGu(stb, enb))
        gptext[st][PTEXT_R_GU] += val;
      else
        gptext[st][PTEXT_R_WC] += val;

      // (   )3<   > 3'
      gptext[st][PTEXT_R] +=
          base01 * Boltzmann(gem.dangle3[en1b][enb][stb]) * gptext[en + 1][PTEXT_R];
      // 5(   )<   > 5'
      gptext[st][PTEXT_R] +=
          base10 * Boltzmann(gem.dangle5[enb][stb][st1b]) * gptext[en + 1][PTEXT_R];
      // .(   ).<   > Terminal mismatch
      gptext[st][PTEXT_R] +=
          base11 * Boltzmann(gem.terminal[en1b][enb][stb][st1b]) * gptext[en + 1][PTEXT_R];
      // .(   ).<(   ) > Left coax
      val = base11 * Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b));
      gptext[st][PTEXT_R] += val * gptext[en + 1][PTEXT_R_GU];
      gptext[st][PTEXT_R] += val * gptext[en + 1][PTEXT_R_WC];

      // (   )<.(   ). > Right coax forward
      gptext[st][PTEXT_R] += base00 * gptext[en + 1][PTEXT_R_RCOAX];
      // (   )<.( * ). > Right coax backward
      gptext[st][PTEXT_R_RCOAX] +=
          base11 * Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b)) * gptext[en + 1][PTEXT_R];

      // (   )(<   ) > Flush coax
      gptext[st][PTEXT_R] +=
          base01 * Boltzmann(gem.stack[en1b][enb][enb ^ 3][stb]) * gptext[en][PTEXT_R_WC];
      if (enb == G || enb == U)
        gptext[st][PTEXT_R] +=
            base01 * Boltzmann(gem.stack[en1b][enb][enb ^ 1][stb]) * gptext[en][PTEXT_R_GU];
    }
  }

  gptext[0][PTEXT_L] = 1;
  for (int en = 1; en < N; ++en) {
    // Case: No pair ending here
    gptext[en][PTEXT_L] += gptext[en - 1][PTEXT_L];
    for (int st = 0; st < en - HAIRPIN_MIN_SZ; ++st) {
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const penergy_t base00 = gpt[st][en][PT_P] * Boltzmann(gem.AuGuPenalty(stb, enb));
      const penergy_t base01 = gpt[st][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(stb, en1b));
      const penergy_t base10 = gpt[st + 1][en][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, enb));
      const penergy_t base11 = gpt[st + 1][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, en1b));
      penergy_t ptextl{1};
      penergy_t ptextlgu{0};
      penergy_t ptextlwc{0};
      penergy_t ptextllcoaxx{0};

      if (st) {
        ptextl = gptext[st - 1][PTEXT_L];
        ptextlgu = gptext[st - 1][PTEXT_L_GU];
        ptextlwc = gptext[st - 1][PTEXT_L_WC];
        ptextllcoaxx = gptext[st - 1][PTEXT_L_LCOAX];
      }

      // <   >(   )
      penergy_t val = base00 * ptextl;
      gptext[en][PTEXT_L] += val;
      if (IsGu(stb, enb))
        gptext[en][PTEXT_L_GU] += val;
      else
        gptext[en][PTEXT_L_WC] += val;

      // <   >(   )3 3'
      gptext[en][PTEXT_L] += base01 * Boltzmann(gem.dangle3[en1b][enb][stb]) * ptextl;
      // <   >5(   ) 5'
      gptext[en][PTEXT_L] += base10 * Boltzmann(gem.dangle5[enb][stb][st1b]) * ptextl;
      // <   >.(   ). Terminal mismatch
      gptext[en][PTEXT_L] += base11 * Boltzmann(gem.terminal[en1b][enb][stb][st1b]) * ptextl;
      // <  (   )>.(   ). Right coax
      val = base11 * Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b));
      gptext[en][PTEXT_L] += val * ptextlgu;
      gptext[en][PTEXT_L] += val * ptextlwc;

      // <  .(   ).>(   ) Left coax forward
      gptext[en][PTEXT_L] += base00 * ptextllcoaxx;
      // <  .( * ).>(   ) Left coax backward
      gptext[en][PTEXT_L_LCOAX] +=
          base11 * Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b)) * ptextl;

      // < (   >)(   ) Flush coax
      gptext[en][PTEXT_L] +=
          base10 * Boltzmann(gem.stack[stb][st1b][enb][stb ^ 3]) * gptext[st][PTEXT_L_WC];
      if (stb == G || stb == U)
        gptext[en][PTEXT_L] +=
            base10 * Boltzmann(gem.stack[stb][st1b][enb][stb ^ 1]) * gptext[st][PTEXT_L_GU];
    }
  }

  assert(fabs(gptext[N - 1][PTEXT_L] - gptext[0][PTEXT_R]) < EP);
}

}  // namespace internal
}  // namespace partition
}  // namespace mrna
