// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#include <energy/energy_model.h>

#include "partition/partition.h"

namespace memerna {
namespace partition {

probabilities_t ComputeProbabilities(const partition_t& partition) {
  const int N = int(partition.p.Size());
  probabilities_t probabilities(std::size_t(N), 0);
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j)
      probabilities[i][j][0] = partition.p[i][j][0] * partition.p[j][i][0] / partition.q;
  return probabilities;
}

}  // namespace partition
}  // namespace memerna
