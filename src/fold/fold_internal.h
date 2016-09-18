#ifndef MEMERNA_FOLD_INTERNAL_H
#define MEMERNA_FOLD_INTERNAL_H

#include "common.h"
#include "base.h"
#include "fold/globals.h"

namespace memerna {
namespace fold {
namespace internal {

// MFE folding related

inline bool ViableFoldingPair(int st, int en) {
  return CanPair(gr[st], gr[en]) &&
      ((en - st - 3 >= HAIRPIN_MIN_SZ && CanPair(gr[st + 1], gr[en - 1])) ||
          (st > 0 && en < int(gr.size() - 1) && CanPair(gr[st - 1], gr[en + 1])));
}

struct cand_t {
  energy_t energy;
  int idx;
};

void ComputeTables0();
void ComputeTables1();
void ComputeTables2();
void ComputeTables3();
void ComputeExterior();
void Traceback();

// Suboptimal folding related:

// Use int16_t here to save memory.
struct index_t {
  int16_t st, en, a;

  index_t() : st(-1), en(-1), a(-1) {}
  index_t(int st_, int en_, int a_) : st(int16_t(st_)), en(int16_t(en_)), a(int16_t(a_)) {
    assert(st_ == st && en_ == en && a == a_);
  }

  bool operator==(const index_t& o) const {
    return st == o.st && en == o.en && a == o.a;
  }

  bool operator!=(const index_t& o) const {
    return !(*this == o);
  }

  bool operator<(const index_t& o) const {
    if (st != o.st) return st < o.st;
    if (en != o.en) return en < o.en;
    return a < o.a;
  }
};

struct ctd_idx_t {
  ctd_idx_t() : idx(-1), ctd(CTD_NA) {}
  ctd_idx_t(int idx_, Ctd ctd_) : idx(int16_t(idx_)), ctd(ctd_) {
    assert(idx_ == idx);
  }

  int16_t idx;
  Ctd ctd;
};

}
}
}

namespace std {
template<>
struct hash<memerna::fold::internal::index_t> {
  size_t operator()(const memerna::fold::internal::index_t& o) const {
    size_t v = 0;
    v ^= hash<int16_t>()(o.st) + 0x9e3779b9 + (v << 6) + (v >> 2);
    v ^= hash<int16_t>()(o.en) + 0x9e3779b9 + (v << 6) + (v >> 2);
    v ^= hash<int16_t>()(o.a) + 0x9e3779b9 + (v << 6) + (v >> 2);
    return v;
  }
};
}

#endif //MEMERNA_FOLD_INTERNAL_H
