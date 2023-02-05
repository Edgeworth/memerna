// Copyright 2016 Eliot Courtney.
#include "structure.h"

#include <cassert>
#include <string>
#include <vector>

namespace mrna {

std::vector<std::string> Structure::Description(int nesting) const {
  std::vector<std::string> desc;
  desc.push_back(fmt::format("{} - {}", nesting, ShortDesc()));
  for (const auto& note : notes_) desc.push_back(" | " + note);
  for (int i = 0; i < static_cast<int>(branches_.size()); ++i)
    desc.push_back(" |-- " + BranchDesc(i));
  for (const auto& branch : branches_) {
    auto branch_desc = branch->Description(nesting + 1);
    desc.insert(desc.end(), branch_desc.begin(), branch_desc.end());
  }
  return desc;
}

std::string HairpinLoopStructure::ShortDesc() const {
  return fmt::format("Hairpin({}, {}) - {}e:{}e", st_, en_, total_energy(), self_energy());
}

std::string TwoLoopStructure::ShortDesc() const {
  return fmt::format(
      "TwoLoop({}, {}, {}, {}) - {}e:{}e", ost_, oen_, ist_, ien_, total_energy(), self_energy());
}

void TwoLoopStructure::AddBranch(std::unique_ptr<Structure> b) {
  assert(branches_.empty());
  Structure::AddBranch(std::move(b));
}

std::string MultiLoopStructure::ShortDesc() const {
  return fmt::format("MultiLoop({}, {}) - {}e:{}e", st_, en_, total_energy(), self_energy());
}

std::string MultiLoopStructure::BranchDesc(int idx) const {
  return fmt::format("{} - {}e - {}", branches_[idx]->ShortDesc(), branch_ctd_[idx].second,
      CtdToName(branch_ctd_[idx].first));
}

std::string StackingStructure::ShortDesc() const {
  return fmt::format("Stacking({}, {}) - {}e:{}e", st_, en_, total_energy(), self_energy());
}

void StackingStructure::AddBranch(std::unique_ptr<Structure> b) {
  assert(branches_.empty());
  Structure::AddBranch(std::move(b));
}

}  // namespace mrna
