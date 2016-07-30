#include "fold.h"
#include "globals.h"

namespace memerna {
namespace fold {

namespace {

array2d_t paired, unpaired;

}

energy_t FoldInternal() {
  int N = int(r.size());
  paired = array2d_t(r.size());
  unpaired = array2d_t(r.size());

  for (int st = 0; st < N; ++st) {
    for (int en = st; en < N; ++en) {
      // Internal loops and bulge loops. TODO: Lyngso's
      //for (int ist = )
    }
  }
}

energy_t Fold(const rna_t& rna, std::unique_ptr<structure::Structure>* s) {
  r = rna;
  return FoldInternal();
}

}
}
