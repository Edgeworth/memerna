#include <stack>
#include "fold/fold.h"
#include "fold/slow_fold.h"
#include "fold/traceback.h"

namespace memerna {
namespace fold {

energy_t Fold() {
  auto arr = ComputeTables_Slow();
  std::stack<std::tuple<int, int, int>> q;
  auto energy = TraceExterior(arr, q);
  TraceStructure(arr, q);
  return energy;
}

}
}
