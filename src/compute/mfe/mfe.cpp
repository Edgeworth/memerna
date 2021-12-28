// Copyright 2016 E.
#include "compute/mfe/mfe.h"

#include <algorithm>
#include <stack>

#include "compute/dp.h"
#include "compute/energy/globals.h"
#include "compute/traceback/traceback.h"

namespace mrna::mfe::internal {

// TODO: get rid of these
using energy::gem;
using mfe::internal::gdp;

void ComputeExterior() {
  const int N = static_cast<int>(gr.size());
  // Exterior loop calculation. There can be no paired base on gext[en].
  gext[N][EXT] = 0;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    gext[st][EXT] = gext[st + 1][EXT];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const auto base00 = gdp[st][en][DP_P] + gem.AuGuPenalty(stb, enb);
      const auto base01 = gdp[st][en - 1][DP_P] + gem.AuGuPenalty(stb, en1b);
      const auto base10 = gdp[st + 1][en][DP_P] + gem.AuGuPenalty(st1b, enb);
      const auto base11 = gdp[st + 1][en - 1][DP_P] + gem.AuGuPenalty(st1b, en1b);
      Energy ext = MAX_E;

      // (   )<   >
      auto val = base00 + gext[en + 1][EXT];
      ext = std::min(ext, val);
      if (IsGu(stb, enb))
        gext[st][EXT_GU] = std::min(gext[st][EXT_GU], val);
      else
        gext[st][EXT_WC] = std::min(gext[st][EXT_WC], val);

      // (   )3<   > 3'
      ext = std::min(ext, base01 + gem.dangle3[en1b][enb][stb] + gext[en + 1][EXT]);
      // 5(   )<   > 5'
      ext = std::min(ext, base10 + gem.dangle5[enb][stb][st1b] + gext[en + 1][EXT]);
      // .(   ).<   > Terminal mismatch
      ext = std::min(ext, base11 + gem.terminal[en1b][enb][stb][st1b] + gext[en + 1][EXT]);
      // .(   ).<(   ) > Left coax
      val = base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b);
      ext = std::min(ext, val + gext[en + 1][EXT_GU]);
      ext = std::min(ext, val + gext[en + 1][EXT_WC]);

      // (   )<.(   ). > Right coax forward
      ext = std::min(ext, base00 + gext[en + 1][EXT_RCOAX]);
      // (   )<.( * ). > Right coax backward
      gext[st][EXT_RCOAX] = std::min(gext[st][EXT_RCOAX],
          base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b) + gext[en + 1][EXT]);

      // (   )(<   ) > Flush coax
      ext = std::min(ext, base01 + gem.stack[en1b][enb][enb ^ 3][stb] + gext[en][EXT_WC]);
      if (enb == G || enb == U)
        ext = std::min(ext, base01 + gem.stack[en1b][enb][enb ^ 1][stb] + gext[en][EXT_GU]);

      gext[st][EXT] = std::min(gext[st][EXT], ext);
    }
  }
}

}  // namespace mrna::mfe::internal
