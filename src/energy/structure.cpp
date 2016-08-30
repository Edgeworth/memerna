#include "structure.h"

namespace memerna {
namespace structure {

const char* CtdToName(Ctd ctd) {
  switch (ctd) {
    case CTD_NA: return "n/a";
    case CTD_UNUSED: return "unused";
    case CTD_3_DANGLE: return "3' dangle";
    case CTD_5_DANGLE: return "5' dangle";
    case CTD_TERMINAL_MISMATCH: return "terminal mismatch";
    case CTD_LEFT_MISMATCH_COAX_WITH_NEXT: return "left mismatch coax with next";
    case CTD_LEFT_MISMATCH_COAX_WITH_PREV: return "left mismatch coax with prev";
    case CTD_RIGHT_MISMATCH_COAX_WITH_NEXT: return "right mismatch coax with next";
    case CTD_RIGHT_MISMATCH_COAX_WITH_PREV: return "right mismatch coax with prev";
    case CTD_FLUSH_COAX_WITH_NEXT: return "flush coax with next";
    case CTD_FLUSH_COAX_WITH_PREV: return "flush coax with prev";
    default: verify_expr(false, "bug");
  }
}

std::vector<std::string> Structure::Description(int nesting) {
  std::vector<std::string> desc;
  desc.push_back(sfmt("%d - %s", nesting, ShortDesc().c_str()));
  for (const auto& note : notes)
    desc.push_back(" | " + note);
  for (int i = 0; i < int(branches.size()); ++i)
    desc.push_back(" |-- " + BranchDesc(i));
  for (const auto& branch : branches) {
    auto branch_desc = branch->Description(nesting + 1);
    desc.insert(desc.end(), branch_desc.begin(), branch_desc.end());
  }
  return desc;
}

void Structure::AddNote(std::string note, ...) {
  va_list l;
  va_start(l, note);
  notes.push_back(vsfmt(note.c_str(), l));
  va_end(l);
}

std::string HairpinLoop::ShortDesc() {
  return sfmt("HairpinLoop(%d, %d) - %de:%de", st, en, GetTotalEnergy(), GetSelfEnergy());
}

std::string InternalLoop::ShortDesc() {
  return sfmt(
      "InternalLoop(%d, %d, %d, %d) - %de:%de",
      ost, oen, ist, ien, GetTotalEnergy(), GetSelfEnergy());
}

std::string MultiLoop::ShortDesc() {
  return sfmt("MultiLoop(%d, %d) - %de:%de", st, en, GetTotalEnergy(), GetSelfEnergy());
}

std::string Stacking::ShortDesc() {
  return sfmt("Stacking(%d, %d) - %de:%de", st, en, GetTotalEnergy(), GetSelfEnergy());
}

}
}
