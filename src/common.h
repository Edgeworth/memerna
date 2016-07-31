#ifndef MEMERNA_COMMON_H
#define MEMERNA_COMMON_H

#include <cassert>
#include <cstdint>
#include <string>
#include <cstdarg>

// Note that changing these might break tests.
#define COMPUTE_CTDS 0
#define USE_HACK_MODEL 1

std::string sfmt(const char* fmt, ...);
std::string vsfmt(const char* fmt, va_list l);

#endif //MEMERNA_COMMON_H
