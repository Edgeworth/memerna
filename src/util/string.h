// Copyright 2021 Eliot Courtney.
#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include "common.h"

namespace mrna {

std::string sgetline(FILE* fp);
std::string sfmt(const char* fmt, ...);
std::string vsfmt(const char* fmt, va_list l);

uint32_t Crc32(const std::string& data);

}  // namespace mrna

#endif  // UTIL_STRING_H_
