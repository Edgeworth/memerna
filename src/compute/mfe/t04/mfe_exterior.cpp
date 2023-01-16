// Copyright 2016 Eliot Courtney.
#include <algorithm>
#include <memory>

#include "compute/dp.h"
#include "compute/energy/t04/model.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::mfe::t04 {

ExtArray MfeExterior(const Primary& r, const erg::t04::Model::Ptr& em, const DpArray& dp) {
  const int N = static_cast<int>(r.size());
  auto ext = ExtArray(r.size() + 1, MAX_E);

  // Exterior loop calculation. There can be no paired base on ext[en].
  ext[N][EXT] = ZERO_E;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    ext[st][EXT] = ext[st + 1][EXT];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const auto base00 = dp[st][en][DP_P] + em->AuGuPenalty(stb, enb);
      const auto base01 = dp[st][en - 1][DP_P] + em->AuGuPenalty(stb, en1b);
      const auto base10 = dp[st + 1][en][DP_P] + em->AuGuPenalty(st1b, enb);
      const auto base11 = dp[st + 1][en - 1][DP_P] + em->AuGuPenalty(st1b, en1b);
      Energy e = MAX_E;

      // (   )<   >
      auto val = base00 + ext[en + 1][EXT];
      e = std::min(e, val);
      if (IsGuPair(stb, enb))
        ext[st][EXT_GU] = std::min(ext[st][EXT_GU], val);
      else
        ext[st][EXT_WC] = std::min(ext[st][EXT_WC], val);

      // (   )3<   > 3'
      e = std::min(e, base01 + em->dangle3[en1b][enb][stb] + ext[en + 1][EXT]);
      // 5(   )<   > 5'
      e = std::min(e, base10 + em->dangle5[enb][stb][st1b] + ext[en + 1][EXT]);
      // .(   ).<   > Terminal mismatch
      e = std::min(e, base11 + em->terminal[en1b][enb][stb][st1b] + ext[en + 1][EXT]);
      // .(   ).<(   ) > Left coax
      val = base11 + em->MismatchCoaxial(en1b, enb, stb, st1b);
      e = std::min(e, val + ext[en + 1][EXT_GU]);
      e = std::min(e, val + ext[en + 1][EXT_WC]);

      // (   )<.(   ). > Right coax forward
      e = std::min(e, base00 + ext[en + 1][EXT_RC]);
      // (   )<.( * ). > Right coax backward
      ext[st][EXT_RC] = std::min(
          ext[st][EXT_RC], base11 + em->MismatchCoaxial(en1b, enb, stb, st1b) + ext[en + 1][EXT]);

      // (   )(<   ) > Flush coax
      e = std::min(e, base01 + em->stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC]);
      if (IsGu(enb))
        e = std::min(e, base01 + em->stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU]);

      ext[st][EXT] = std::min(ext[st][EXT], e);
    }
  }

  return ext;
}

}  // namespace mrna::mfe::t04
