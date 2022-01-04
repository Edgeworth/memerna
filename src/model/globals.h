// Copyright 2016 Eliot Courtney.
#ifndef MODEL_GLOBALS_H_
#define MODEL_GLOBALS_H_

#include "model/model.h"
#include "model/primary.h"

namespace mrna {

extern Primary gr;
void SetGlobalState(const Primary& r);

}  // namespace mrna

#endif  // MODEL_GLOBALS_H_
