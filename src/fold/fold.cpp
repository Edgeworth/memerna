#include <stack>
#include "fold/fold.h"
#include "fold/slow_fold.h"
#include "fold/traceback.h"
#include "fold/fold1.h"

namespace memerna {
namespace fold {

template<typename T>
energy_t FoldInternal(T ComputeTables) {
  auto arr = ComputeTables();
  std::stack<std::tuple<int, int, int>> q;
  auto energy = TraceExterior(arr, q);
  TraceStructure(arr, q);
  return energy;
}

energy_t Fold() {
  return FoldInternal(ComputeTablesSlow);
}

energy_t Fold1() {
  return FoldInternal(ComputeTables1);
}

}
}
