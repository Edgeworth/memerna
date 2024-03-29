// Copyright 2022 E.
#ifndef MODELS_T04_ENERGY_PRECOMP_H_
#define MODELS_T04_ENERGY_PRECOMP_H_
#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "model/base.h"
#include "model/energy.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "util/error.h"

namespace mrna::md::t04 {

inline constexpr int MAX_SPECIAL_HAIRPIN_SZ = 6;

// This is templated because the partition function wants to use it with a different type.
template <typename T>
struct HairpinPrecomp {
  explicit HairpinPrecomp(T init) {
    std::fill(special, special + sizeof(special) / sizeof(special[0]), init);
  }

  T special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c{0};
};

template <typename HairpinPrecomp, bool is_boltz, typename EM>
std::vector<HairpinPrecomp> PrecomputeHairpin(const Primary& r, const EM& em, auto init) {
  std::vector<HairpinPrecomp> pc(r.size(), HairpinPrecomp(init));
  std::string rna_str = r.ToSeq();
  for (const auto& hairpinpair : em.hairpin) {
    const auto& str = hairpinpair.first;
    verify(str.size() - 2 <= MAX_SPECIAL_HAIRPIN_SZ, "need to increase MAX_SPECIAL_HAIRPIN_SZ");
    auto energy = hairpinpair.second;
    auto augu = em.AuGuPenalty(CharToBase(*str.begin()).value(), CharToBase(*str.rbegin()).value());
    if constexpr (is_boltz)
      energy *= augu;
    else
      energy += augu;

    auto pos = rna_str.find(str, 0);
    while (pos != std::string::npos) {
      pc[pos].special[str.size() - 2] = energy;
      pos = rna_str.find(str, pos + 1);
    }
  }
  const int N = static_cast<int>(r.size());
  pc[N - 1].num_c = static_cast<int>(r[N - 1] == C);
  for (int i = N - 2; i >= 0; --i)
    if (r[i] == C) pc[i].num_c = pc[i + 1].num_c + 1;
  return pc;
}

class Precomp {
 public:
  Energy augubranch[4][4]{};
  Energy min_mismatch_coax{};
  Energy min_flush_coax{};
  Energy min_twoloop_not_stack{};

  std::vector<HairpinPrecomp<Energy>> hairpin;

  Precomp(Primary r, Model::Ptr em);

  [[nodiscard]] const Model& em() const { return *em_; }

  [[nodiscard]] Energy TwoLoop(int ost, int oen, int ist, int ien) const;
  [[nodiscard]] Energy Hairpin(int st, int en) const;

 private:
  Primary r_;
  Model::Ptr em_;

  void PrecomputeData();
};

}  // namespace mrna::md::t04

#endif  // MODELS_T04_ENERGY_PRECOMP_H_
