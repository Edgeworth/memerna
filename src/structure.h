#ifndef MEMERNA_STRUCTURE_H
#define MEMERNA_STRUCTURE_H

#include <memory>
#include "common.h"
#include "energy.h"

namespace memerna {
namespace structure {

class Structure {
public:
  Structure() = default;
  Structure(const Structure&) = delete;

  Structure& operator=(const Structure&) = delete;

  void AddNote(const std::string& note) {
    notes.push_back(note);
  }
  std::vector<std::string> Description(int nesting);
  virtual std::string ShortDesc() = 0;
  virtual std::vector<std::string> BodyDesc() = 0;

  void SetEnergy(energy::energy_t e) {energy = e;}
private:
  energy::energy_t energy;
  std::vector<std::string> notes;
};


class HairpinLoop : public Structure {
public:
  HairpinLoop(int st, int en) : st(st), en(en) {}

  std::string ShortDesc();
  std::vector<std::string> BodyDesc();
private:
  int st, en;
};

class BulgeLoop : public Structure {
public:
  BulgeLoop(int ost, int oen, int ist, int ien) : ost(ost), oen(oen), ist(ist), ien(ien) {}

  void SetBranch(std::unique_ptr<Structure> c) {
    branch = std::move(c);
  }
  std::string ShortDesc();
  std::vector<std::string> BodyDesc();
private:
  int ost, oen, ist, ien;
  std::unique_ptr<Structure> branch;
};

class InternalLoop : public Structure {
public:
  InternalLoop(int ost, int oen, int ist, int ien) : ost(ost), oen(oen), ist(ist), ien(ien) {}

  void SetBranch(std::unique_ptr<Structure> b) {
    branch = std::move(b);
  }
  std::string ShortDesc();
  std::vector<std::string> BodyDesc();
private:
  int ost, oen, ist, ien;
  std::unique_ptr<Structure> branch;
};

enum class CtdType {
  UNUSED,
  STOLEN,  // Used by another branch.
  NOEXIST,
  LEFT_DANGLE,
  RIGHT_DANGLE,
  TERMINAL_MISMATCH,
  COAXIAL_MIDDLE,  // Involved in a mismatch
  COAXIAL_END
};

struct ctd_t {

  // Also store if this branch is involved in a flush coaxial stack.
};

class MultiLoop : public Structure {
public:
  MultiLoop(int st, int en) : st(st), en(en) {}

  void AddBranch(std::unique_ptr<Structure> b) {
    branches.push_back(std::move(b));
  }
  std::string ShortDesc();
  std::vector<std::string> BodyDesc();
private:
  int st, en;
  std::vector<std::unique_ptr<Structure>> branches;
  // TODO ctds.
};

class Stacking : public Structure {
public:
  Stacking(int st, int en) : st(st), en(en) {}

  void SetBranch(std::unique_ptr<Structure> b) {
    branch = std::move(b);
  }
  std::string ShortDesc();
  std::vector<std::string> BodyDesc();
private:
  int st, en;
  std::unique_ptr<Structure> branch;
};

}
}

#endif //MEMERNA_PARSING_H
