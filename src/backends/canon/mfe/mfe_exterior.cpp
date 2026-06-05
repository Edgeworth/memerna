// Copyright 2026 Eliot Courtney.
#include "backends/canon/mfe/mfe_exterior.h"

#include <algorithm>

#include "api/energy/energy_cfg.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
Energy MfeExterior(
    const Primary& r, const Model::Ptr& m, DpState& state, const PseudofreeCfgT& pf) {
  const int N = r.size();
  pf.Verify(r);

  state.ext = ExtArray(r.size() + 1, MAX_E);
  auto& [dp, ext] = state;

  ext[N][EXT] = ZERO_E;
  for (int st = N - 1; st >= 0; --st) {
    ext[st][EXT] = ext[st + 1][EXT] + pf.Unpaired(st);
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const auto base00 = dp[st][en][DP_P] + m->AuGuPenalty(stb, enb);
      const auto base01 = dp[st][en - 1][DP_P] + m->AuGuPenalty(stb, en1b);
      const auto base10 = dp[st + 1][en][DP_P] + m->AuGuPenalty(st1b, enb);
      const auto base11 = dp[st + 1][en - 1][DP_P] + m->AuGuPenalty(st1b, en1b);

      auto val = base00 + ext[en + 1][EXT];

      if constexpr (CFG.UseD2()) {
        if (st != 0 && en != N - 1)
          val += m->terminal[enb][r[en + 1]][r[st - 1]][stb];
        else if (en != N - 1)
          val += m->dangle3[enb][r[en + 1]][stb];
        else if (st != 0)
          val += m->dangle5[enb][r[st - 1]][stb];
      }

      Energy e = val;

      if constexpr (CFG.UseDangleMismatch()) {
        e = std::min(e, base01 + m->dangle3[en1b][enb][stb] + pf.Unpaired(en) + ext[en + 1][EXT]);
        e = std::min(e, base10 + m->dangle5[enb][stb][st1b] + pf.Unpaired(st) + ext[en + 1][EXT]);
        e = std::min(e,
            base11 + m->terminal[en1b][enb][stb][st1b] + pf.Unpaired(st) + pf.Unpaired(en) +
                ext[en + 1][EXT]);
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // .(   ).<(   ) > Left coax
        for (int piv = en + HAIRPIN_MIN_SZ + 2; piv < N; ++piv) {
          const Base prb = r[en + 1];
          const Base pb = r[piv];
          e = std::min(e,
              base11 + dp[en + 1][piv][DP_P] + m->AuGuPenalty(prb, pb) +
                  m->MismatchCoaxial(en1b, enb, stb, st1b) + pf.Unpaired(st) + pf.Unpaired(en) +
                  ext[piv + 1][EXT]);
        }

        // (   )<.(   ). > Right coax forward
        for (int piv = en + HAIRPIN_MIN_SZ + 4; piv < N; ++piv) {
          const Base prb = r[en + 1];
          const Base pr1b = r[en + 2];
          const Base pl1b = r[piv - 1];
          const Base pb = r[piv];
          e = std::min(e,
              base00 + dp[en + 2][piv - 1][DP_P] + m->AuGuPenalty(pr1b, pl1b) +
                  m->MismatchCoaxial(pl1b, pb, prb, pr1b) + pf.Unpaired(en + 1) + pf.Unpaired(piv) +
                  ext[piv + 1][EXT]);
        }

        // (   )(<   ) > Flush coax
        for (int piv = en + HAIRPIN_MIN_SZ + 1; piv < N; ++piv) {
          const Base pb = r[piv];
          e = std::min(e,
              base01 + dp[en][piv][DP_P] + m->AuGuPenalty(enb, pb) + m->stack[en1b][enb][pb][stb] +
                  ext[piv + 1][EXT]);
        }
      }

      ext[st][EXT] = std::min(ext[st][EXT], e);
    }
  }

  return ext[0][EXT];
}

#define MRNA_INSTANTIATE(CFG, PF) \
  template Energy MfeExterior<CFG, PF>(const Primary&, const Model::Ptr&, DpState&, const PF&);
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::canon
