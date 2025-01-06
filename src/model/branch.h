// Copyright 2016 Eliot Courtney.
#ifndef MODEL_BRANCH_H_
#define MODEL_BRANCH_H_

#include <deque>
#include <vector>

#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna {

// For each base in `s`, the value at the index of that base will contain the
// number of child branches in the loop that that base is a child of. For
// example, if the base is part of a base pair inside a multiloop with 2 inner
// branches, the value will be 2. N.B. that bases in the exterior loop will
// always have their value set to 2,  the exterior loop is treated as a
// multiloop.
// N.B. The value at the index of the right side of a branch will contain
// the number of children that branch has. The value at the index of the left
// side of a branch will contain the number of siblings that branch has.
std::vector<int> GetBranchCounts(const Secondary& s);

// Takes the branch representation of ctds in `branch_ctd` for `branches` branches and
// writes it in per-base representation to `ctds`.
void AddBranchCtdsToBaseCtds(
    const std::deque<int>& branches, const BranchCtd& branch_ctd, Ctds* ctd);

int MaxNumContiguous(const Primary& r);

}  // namespace mrna

#endif  // MODEL_BRANCH_H_
