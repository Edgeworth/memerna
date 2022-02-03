// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_STRUCTURE_H_
#define COMPUTE_ENERGY_STRUCTURE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/branch.h"
#include "model/ctd.h"
#include "model/model.h"
#include "util/error.h"

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

  virtual void AddBranch(std::unique_ptr<Structure> b) { branches_.push_back(std::move(b)); }
  virtual std::string BranchDesc(int idx) const { return branches_[idx]->ShortDesc(); }

  void set_self_energy(Energy e) { self_energy_ = e; }
  void set_total_energy(Energy e) { total_energy_ = e; }
  Energy self_energy() const { return self_energy_; }
  Energy total_energy() const { return total_energy_; }

 protected:
  std::vector<std::unique_ptr<Structure>> branches_;

 private:
  Energy self_energy_;
  Energy total_energy_;
  std::vector<std::string> notes_;
};

class HairpinLoopStructure : public Structure {
 public:
  HairpinLoopStructure(int st, int en) : st_(st), en_(en) {}

  void AddBranch(std::unique_ptr<Structure>) override { error("invalid operation"); }
  std::string ShortDesc() const override;

 private:
  int st_, en_;
};

class TwoLoopStructure : public Structure {
 public:
  TwoLoopStructure(int ost, int oen, int ist, int ien)
      : ost_(ost), oen_(oen), ist_(ist), ien_(ien) {}

  void AddBranch(std::unique_ptr<Structure> b) override;
  std::string ShortDesc() const override;

 private:
  int ost_, oen_, ist_, ien_;
};

class MultiLoopStructure : public Structure {
 public:
  MultiLoopStructure(int st, int en) : st_(st), en_(en) {}

  void AddCtd(Ctd ctd, Energy ctd_energy) { branch_ctd_.emplace_back(ctd, ctd_energy); }
  std::string BranchDesc(int idx) const override;
  std::string ShortDesc() const override;

 private:
  int st_, en_;
  BranchCtd branch_ctd_;
};

class StackingStructure : public Structure {
 public:
  StackingStructure(int st, int en) : st_(st), en_(en) {}

  void AddBranch(std::unique_ptr<Structure> b) override;
  std::string ShortDesc() const override;

 private:
  int st_, en_;
};

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_STRUCTURE_H_
