// Copyright 2016 E.
#ifndef FOLD_FOLD_H_
#define FOLD_FOLD_H_

#include "common.h"
#include "mfe/fold_globals.h"
#include "model/base.h"
#include "model/globals.h"

namespace mrna::fold {

namespace internal {

// MFE folding related

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

}  // namespace internal

typedef std::function<void(const computed_t&)> SuboptimalCallback;

}  // namespace mrna::fold

#endif  // FOLD_FOLD_H_
