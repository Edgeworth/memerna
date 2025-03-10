// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna::md::base {

Energy MfeExterior(const Primary& r, const Model::Ptr& m, DpState& state) {
  const int N = static_cast<int>(r.size());

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::D2,
          erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m->cfg());
  m->pf.Verify(r);

  state.ext = ExtArray(r.size() + 1, MAX_E);
  auto& [dp, ext] = state;

  // Exterior loop calculation. There can be no paired base on ext[en].
  ext[N][EXT] = ZERO_E;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    ext[st][EXT] = ext[st + 1][EXT] + m->pf.Unpaired(st);
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const auto base00 = dp[st][en][DP_P] + m->AuGuPenalty(stb, enb);
      const auto base01 = dp[st][en - 1][DP_P] + m->AuGuPenalty(stb, en1b);
      const auto base10 = dp[st + 1][en][DP_P] + m->AuGuPenalty(st1b, enb);
      const auto base11 = dp[st + 1][en - 1][DP_P] + m->AuGuPenalty(st1b, en1b);
      Energy e = MAX_E;

      // (   )<   >
      auto val = base00 + ext[en + 1][EXT];

      if (m->cfg().UseD2()) {
        if (st != 0 && en != N - 1) {
          // (   )<   > Terminal mismatch - U
          val += m->terminal[enb][r[en + 1]][r[st - 1]][stb];
        } else if (en != N - 1) {
          // (   )<3   > 3' - U
          val += m->dangle3[enb][r[en + 1]][stb];
        } else if (st != 0) {
          // 5(   )<   > 5' - U
          val += m->dangle5[enb][r[st - 1]][stb];
        }
      }

      e = std::min(e, val);
      if (IsGuPair(stb, enb))
        ext[st][EXT_GU] = std::min(ext[st][EXT_GU], val);
      else
        ext[st][EXT_WC] = std::min(ext[st][EXT_WC], val);

      if (m->cfg().UseDangleMismatch()) {
        // (   )3<   > 3'
        e = std::min(
            e, base01 + m->dangle3[en1b][enb][stb] + m->pf.Unpaired(en) + ext[en + 1][EXT]);
        // 5(   )<   > 5'
        e = std::min(
            e, base10 + m->dangle5[enb][stb][st1b] + m->pf.Unpaired(st) + ext[en + 1][EXT]);
        // .(   ).<   > Terminal mismatch
        e = std::min(e,
            base11 + m->terminal[en1b][enb][stb][st1b] + m->pf.Unpaired(st) + m->pf.Unpaired(en) +
                ext[en + 1][EXT]);
      }

      if (m->cfg().UseCoaxialStacking()) {
        // .(   ).<(   ) > Left coax
        val = base11 + m->MismatchCoaxial(en1b, enb, stb, st1b) + m->pf.Unpaired(st) +
            m->pf.Unpaired(en);
        e = std::min(e, val + ext[en + 1][EXT_GU]);
        e = std::min(e, val + ext[en + 1][EXT_WC]);

        // (   )<.(   ). > Right coax forward
        e = std::min(e, base00 + ext[en + 1][EXT_RC]);
        // (   )<.( * ). > Right coax backward
        ext[st][EXT_RC] = std::min(ext[st][EXT_RC],
            base11 + m->MismatchCoaxial(en1b, enb, stb, st1b) + m->pf.Unpaired(st) +
                m->pf.Unpaired(en) + ext[en + 1][EXT]);

        // (   )(<   ) > Flush coax
        e = std::min(e, base01 + m->stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC]);
        if (IsGu(enb))
          e = std::min(e, base01 + m->stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU]);
      }

      ext[st][EXT] = std::min(ext[st][EXT], e);
    }
  }

  return ext[0][EXT];
}
}  // namespace mrna::md::base
