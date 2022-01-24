// Copyright 2016 Eliot Courtney.
#include "structure.h"

#include <string>
#include <vector>

#include "util/string.h"

namespace mrna::energy {

const char* CtdToName(Ctd ctd) {
  switch (ctd) {
  case CTD_NA: return "n/a";
  case CTD_UNUSED: return "unused";
  case CTD_3_DANGLE: return "3' dangle";
  case CTD_5_DANGLE: return "5' dangle";
  case CTD_MISMATCH: return "terminal mismatch";
  case CTD_LCOAX_WITH_NEXT: return "left mismatch coax with next";
  case CTD_LCOAX_WITH_PREV: return "left mismatch coax with prev";
  case CTD_RCOAX_WITH_NEXT: return "right mismatch coax with next";
  case CTD_RCOAX_WITH_PREV: return "right mismatch coax with prev";
  case CTD_FCOAX_WITH_NEXT: return "flush coax with next";
  case CTD_FCOAX_WITH_PREV: return "flush coax with prev";
  default: bug();
  }
}

std::vector<std::string> Structure::Description(int nesting) const {
  std::vector<std::string> desc;
  desc.push_back(sfmt("%d - %s", nesting, ShortDesc().c_str()));
  for (const auto& note : notes_) desc.push_back(" | " + note);
  for (int i = 0; i < static_cast<int>(branches_.size()); ++i)
    desc.push_back(" |-- " + BranchDesc(i));
  for (const auto& branch : branches_) {
    auto branch_desc = branch->Description(nesting + 1);
    desc.insert(desc.end(), branch_desc.begin(), branch_desc.end());
  }
  return desc;
}

void Structure::AddNote(std::string note, ...) {
  va_list l;
  va_start(l, note);
  notes_.push_back(vsfmt(note.c_str(), l));
  va_end(l);
}

std::string HairpinLoopStructure::ShortDesc() const {
  return sfmt("Hairpin(%d, %d) - %de:%de", st_, en_, total_energy(), self_energy());
}

std::string TwoLoopStructure::ShortDesc() const {
  return sfmt(
      "TwoLoop(%d, %d, %d, %d) - %de:%de", ost_, oen_, ist_, ien_, total_energy(), self_energy());
}

void TwoLoopStructure::AddBranch(std::unique_ptr<Structure> b) {
  assert(branches_.empty());
  Structure::AddBranch(std::move(b));
}

std::string MultiLoopStructure::ShortDesc() const {
  return sfmt("MultiLoop(%d, %d) - %de:%de", st_, en_, total_energy(), self_energy());
}

std::string MultiLoopStructure::BranchDesc(int idx) const {
  return sfmt("%s - %de - %s", branches_[idx]->ShortDesc().c_str(), branch_ctd_[idx].second,
      CtdToName(branch_ctd_[idx].first));
}

std::string StackingStructure::ShortDesc() const {
  return sfmt("Stacking(%d, %d) - %de:%de", st_, en_, total_energy(), self_energy());
}

void StackingStructure::AddBranch(std::unique_ptr<Structure> b) {
  assert(branches_.empty());
  Structure::AddBranch(std::move(b));
}

}  // namespace mrna::energy
