#include "common.h"

namespace {
const int BUF_SIZE = 1024;
}

std::string sfmt(const char* fmt, ...) {
  va_list l;
  va_start(l, fmt);
  std::string res = vsfmt(fmt, l);
  va_end(l);
  return res;
}

std::string vsfmt(const char* fmt, va_list l) {
  char buffer[BUF_SIZE];
  int res = vsnprintf(buffer, BUF_SIZE, fmt, l);
  assert(res >= 0 && res < BUF_SIZE);
  return buffer;
}
