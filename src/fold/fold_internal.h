#ifndef MEMERNA_FOLD_INTERNAL_H
#define MEMERNA_FOLD_INTERNAL_H

#include "common.h"
#include "base.h"
#include "constants.h"
#include "fold/globals.h"

namespace memerna {
namespace fold {
namespace internal {

inline bool ViableFoldingPair(int st, int en) {
  return CanPair(gr[st], gr[en]) &&
      ((en - st - 3 >= constants::HAIRPIN_MIN_SZ && CanPair(gr[st + 1], gr[en - 1])) ||
          (st > 0 && en < int(gr.size() - 1) && CanPair(gr[st - 1], gr[en + 1])));
}

struct cand_t {
  energy_t energy;
  int idx;
};

struct index_t {
  int st, en, a;

  index_t() = default;
  index_t(int st_, int en_, int a_) : st(st_), en(en_), a(a_) {}

  bool operator==(const index_t& o) const {
    return st == o.st && en == o.en && a == o.a;
  }

  bool operator!=(const index_t& o) const {
    return !(*this == o);
  }

  bool operator<(const index_t& o) const {
    if (st != o.st) return st < o.st;
    if (en != o.en) return en < o.en;
    if (a != o.a) return a < o.a;
    return false;
  }
};

void ComputeTables0();
void ComputeTables1();
void ComputeTables2();
void ComputeTables3();
void ComputeExterior();
void Traceback();

}
}
}

namespace std {
template<>
struct hash<memerna::fold::internal::index_t> {
  size_t operator()(const memerna::fold::internal::index_t& o) const {
    size_t v = 0;
    v ^= hash<int>()(o.st) + 0x9e3779b9 + (v << 6) + (v >> 2);
    v ^= hash<int>()(o.en) + 0x9e3779b9 + (v << 6) + (v >> 2);
    v ^= hash<int>()(o.a) + 0x9e3779b9 + (v << 6) + (v >> 2);
    return v;
  }
};
}

#endif //MEMERNA_FOLD_INTERNAL_H
