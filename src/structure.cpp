#include "structure.h"

namespace memerna {
namespace structure {

std::vector<std::string> Structure::Description(int nesting) {
  std::vector<std::string> desc;
  desc.push_back(sfmt("%d - %de - %s", nesting, energy, ShortDesc().c_str()));
  for (const auto& note : notes)
    desc.push_back(" | " + note);
  for (const auto& branch : branches)
    desc.push_back(" |-- " + branch->ShortDesc());
  for (const auto& branch : branches) {
    auto branch_desc = branch->Description(nesting + 1);
    desc.insert(desc.end(), branch_desc.begin(), branch_desc.end());
  }
  return desc;
}

void Structure::AddNote(const std::string& note, ...) {
  va_list l;
  va_start(l, note);
  notes.push_back(vsfmt(note.c_str(), l));
  va_end(l);
}

std::string HairpinLoop::ShortDesc() {
  return sfmt("HairpinLoop(%d, %d)", st, en);
}

std::string InternalLoop::ShortDesc() {
  return sfmt("InternalLoop(%d, %d, %d, %d)", ost, oen, ist, ien);
}

std::string MultiLoop::ShortDesc() {
  return sfmt("MultiLoop(%d, %d)", st, en);
}

std::string Stacking::ShortDesc() {
  return sfmt("Stacking(%d, %d)", st, en);
}

}
}
