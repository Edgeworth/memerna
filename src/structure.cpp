#include "structure.h"

namespace memerna {
namespace structure {

std::vector<std::string> Structure::Description(int nesting) {
  std::vector<std::string> desc;
  desc.push_back(sfmt("| %d %s - %d", nesting, ShortDesc().c_str(), energy));
  for (const auto& note : notes)
    desc.push_back("| " + note);
  for (const auto& body_desc : BodyDesc())
    desc.push_back("|- " + body_desc);
  return desc;
}

std::string HairpinLoop::ShortDesc() {
  return sfmt("HairpinLoop(%d, %d)", st, en);
}

std::vector<std::string> HairpinLoop::BodyDesc() {
  return {};
}

std::string BulgeLoop::ShortDesc() {
  return sfmt("BulgeLoop(%d, %d, %d, %d)", ost, oen, ist, ien);
}

std::vector<std::string> BulgeLoop::BodyDesc() {
  return {};
}

std::string InternalLoop::ShortDesc() {
  return sfmt("InternalLoop(%d, %d, %d, %d)", ost, oen, ist, ien);
}

std::vector<std::string> InternalLoop::BodyDesc() {
  return {};
}

std::string MultiLoop::ShortDesc() {
  return sfmt("MultiLoop(%d, %d)", st, en);
}

std::vector<std::string> MultiLoop::BodyDesc() {
  return {};
}

std::string Stacking::ShortDesc() {
  return sfmt("Stacking(%d, %d)", st, en);
}

std::vector<std::string> Stacking::BodyDesc() {
  return {};
}

}
}
