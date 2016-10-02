// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef MEMERNA_FOLD_INTERNAL_H
#define MEMERNA_FOLD_INTERNAL_H

#include "base.h"
#include "common.h"
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

  bool operator==(const index_t& o) const { return st == o.st && en == o.en && a == o.a; }
  bool operator!=(const index_t& o) const { return !(*this == o); }
  bool operator<(const index_t& o) const {
    if (st != o.st) return st < o.st;
    if (en != o.en) return en < o.en;
    return a < o.a;
  }
};

struct ctd_idx_t {
  ctd_idx_t() : idx(-1), ctd(CTD_NA) {}
  ctd_idx_t(int idx_, Ctd ctd_) : idx(int16_t(idx_)), ctd(ctd_) { assert(idx_ == idx); }

  int16_t idx;
  Ctd ctd;
};
}
}
}

#endif  // MEMERNA_FOLD_INTERNAL_H
