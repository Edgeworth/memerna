// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_STRUCTURE_H_
#define COMPUTE_ENERGY_STRUCTURE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/energy/internal.h"
#include "util/macros.h"

namespace mrna::energy {

const char* CtdToName(Ctd ctd);

class Structure {
 public:
  Structure() = default;
  Structure(const Structure&) = delete;
  Structure& operator=(const Structure&) = delete;

  virtual ~Structure() = default;

  // This is not a reference because of varargs.
  void AddNote(std::string note, ...);
  std::vector<std::string> Description(int nesting = 0) const;
  virtual std::string ShortDesc() const = 0;

  virtual void AddBranch(std::unique_ptr<Structure> b) { branches.push_back(std::move(b)); }
  virtual std::string BranchDesc(int idx) const { return branches[idx]->ShortDesc(); }

  void SetSelfEnergy(Energy e) { self_energy = e; }
  void SetTotalEnergy(Energy e) { total_energy = e; }
  Energy GetSelfEnergy() const { return self_energy; }
  Energy GetTotalEnergy() const { return total_energy; }

 protected:
  std::vector<std::unique_ptr<Structure>> branches;

 private:
  Energy self_energy;
  Energy total_energy;
  std::vector<std::string> notes;
};

class HairpinLoopStructure : public Structure {
 public:
  HairpinLoopStructure(int st_, int en_) : st(st_), en(en_) {}

  void AddBranch(std::unique_ptr<Structure>) { error("invalid operation"); }
  std::string ShortDesc() const;

 private:
  int st, en;
};

class InternalLoopStructure : public Structure {
 public:
  InternalLoopStructure(int ost_, int oen_, int ist_, int ien_)
      : ost(ost_), oen(oen_), ist(ist_), ien(ien_) {}

  void AddBranch(std::unique_ptr<Structure> b);
  std::string ShortDesc() const;

 private:
  int ost, oen, ist, ien;
};

class MultiLoopStructure : public Structure {
 public:
  MultiLoopStructure(int st_, int en_) : st(st_), en(en_) {}

  void AddCtd(Ctd ctd, Energy ctd_energy) { branch_ctds.emplace_back(ctd, ctd_energy); }
  std::string BranchDesc(int idx) const;
  std::string ShortDesc() const;

 private:
  int st, en;
  internal::BranchCtd branch_ctds;
};

class StackingStructure : public Structure {
 public:
  StackingStructure(int st_, int en_) : st(st_), en(en_) {}

  void AddBranch(std::unique_ptr<Structure> b);
  std::string ShortDesc() const;

 private:
  int st, en;
};

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_STRUCTURE_H_
