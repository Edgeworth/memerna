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

  virtual ~Structure() = default;

  void AddNote(const std::string& note, ...);

  std::vector<std::string> Description(int nesting = 0);
  virtual std::string ShortDesc() = 0;

  virtual void AddBranch(std::unique_ptr<Structure> b) {
    branches.push_back(std::move(b));
  }

  void SetEnergy(energy::energy_t e) {energy = e;}

  energy::energy_t GetEnergy() {return energy;}

protected:
  std::vector<std::unique_ptr<Structure>> branches;

private:
  energy::energy_t energy;
  std::vector<std::string> notes;
};


class HairpinLoop : public Structure {
public:
  HairpinLoop(int st, int en) : st(st), en(en) {}

  void AddBranch(std::unique_ptr<Structure>) {
    assert(false);
  }

  std::string ShortDesc();
private:
  int st, en;
};

class InternalLoop : public Structure {
public:
  InternalLoop(int ost, int oen, int ist, int ien) : ost(ost), oen(oen), ist(ist), ien(ien) {}

  void AddBranch(std::unique_ptr<Structure> b) {
    assert(branches.empty());
    Structure::AddBranch(std::move(b));
  }

  std::string ShortDesc();
private:
  int ost, oen, ist, ien;
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


  std::string ShortDesc();
private:
  int st, en;
  // TODO ctds.
};

class Stacking : public Structure {
public:
  Stacking(int st, int en) : st(st), en(en) {}

  void AddBranch(std::unique_ptr<Structure> b) {
    assert(branches.empty());
    Structure::AddBranch(std::move(b));
  }

  std::string ShortDesc();
private:
  int st, en;
};

}
}

#endif //MEMERNA_PARSING_H
