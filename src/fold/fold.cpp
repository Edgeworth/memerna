#include <stack>
#include <climits>
#include "parsing.h"
#include "fold/fold.h"
#include "fold/suboptimal0.h"

namespace memerna {
namespace fold {

using namespace constants;
using namespace energy;

constexpr context_options_t::TableAlg context_options_t::TABLE_ALGS[];
constexpr context_options_t::SuboptimalAlg context_options_t::SUBOPTIMAL_ALGS[];

context_options_t ContextOptionsFromArgParse(const ArgParse& argparse) {
  context_options_t options;
  auto opt = argparse.GetOption("alg");
  if (opt == "0")
    options.table_alg = context_options_t::TableAlg::ZERO;
  else if (opt == "1")
    options.table_alg = context_options_t::TableAlg::ONE;
  else if (opt == "2")
    options.table_alg = context_options_t::TableAlg::TWO;
  else if (opt == "3")
    options.table_alg = context_options_t::TableAlg::THREE;
  else
    verify_expr(false, "unknown fold option");
  options.subopt_energy = atoi(argparse.GetOption("subopt-delta").c_str());
  options.subopt_num = atoi(argparse.GetOption("subopt-num").c_str());
  return options;
}

Context::Context(const primary_t& r_, const energy::EnergyModel& em_)
    : r(r_), em(em_), options(), N(int(r_.size())), pc(internal::PrecomputeData(r_, em_)),
      arr(r_.size() + 1), exterior(r_.size()) {}

Context::Context(const primary_t& r_, const energy::EnergyModel& em_, context_options_t options_)
    : r(r_), em(em_), options(options_), N(int(r_.size())),
      pc(internal::PrecomputeData(r_, em_)), arr(r_.size() + 1), exterior(r_.size() + 1) {}

void Context::ComputeTables() {
  switch (options.table_alg) {
    case context_options_t::TableAlg::ZERO:
      ComputeTables0();
      break;
    case context_options_t::TableAlg::ONE:
      ComputeTables1();
      break;
    case context_options_t::TableAlg::TWO:
      ComputeTables2();
      break;
    case context_options_t::TableAlg::THREE:
      ComputeTables3();
      break;
  }
}

computed_t Context::Fold() {
  ComputeTables();
  ComputeExterior();
  return Traceback();
}

std::vector<computed_t> Context::Suboptimal() {
  // TODO potential recomputation.
  ComputeTables();
  ComputeExterior();
  energy_t max_energy = exterior[0][EXT] + options.subopt_energy;
  int max_structures = options.subopt_num;
  if (options.subopt_energy < 0) max_energy = constants::CAP_E;
  if (options.subopt_num < 0) max_structures = INT_MAX;
  switch (options.suboptimal_alg) {
    case context_options_t::SuboptimalAlg::ZERO:
      return Suboptimal0(*this, max_energy, max_structures).Run();
    default:
      verify_expr(false, "bug");
  }
}

energy_t Context::FastTwoLoop(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0)
    return em.stack[r[ost]][r[ist]][r[ien]][r[oen]];
  if (toplen == 0 || botlen == 0)
    return em.Bulge(r, ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return em.internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return em.internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return em.internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return em.internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];

  static_assert(TWOLOOP_MAX_SZ <= EnergyModel::INITIATION_CACHE_SZ, "initiation cache not large enough");
  energy_t energy =
      em.internal_init[toplen + botlen] +
          std::min(std::abs(toplen - botlen) * em.internal_asym, NINIO_MAX_ASYM);

  energy += em.InternalLoopAuGuPenalty(r[ost], r[oen]);
  energy += em.InternalLoopAuGuPenalty(r[ist], r[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += em.internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        em.internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += em.internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        em.internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];

  return energy;
}

energy_t Context::FastHairpin(int st, int en) {
  int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= internal::hairpin_precomp_t::MAX_SPECIAL_HAIRPIN_SZ && pc.hairpin[st].special[length] != MAX_E)
    return pc.hairpin[st].special[length];
  base_t stb = r[st], st1b = r[st + 1], en1b = r[en - 1], enb = r[en];
  energy_t energy = em.HairpinInitiation(length) + em.AuGuPenalty(stb, enb);

  bool all_c = pc.hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c)
      energy += em.hairpin_c3_loop;
    return energy;
  }
  energy += em.terminal[r[st]][st1b][en1b][r[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy += em.hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G)
    energy += em.hairpin_gg_first_mismatch;
  if (all_c)
    energy += em.hairpin_all_c_a * length + em.hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r[st - 1] == G && r[st - 2] == G)
    energy += em.hairpin_special_gu_closure;

  return energy;
}

}
}
