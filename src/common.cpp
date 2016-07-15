#include "common.h"

namespace {
const int BUF_SIZE = 1024;
}

std::string sfmt(const char* fmt, ...) {
  va_list l;
  va_start(l, fmt);
  char buffer[BUF_SIZE];
  int res = vsnprintf(buffer, BUF_SIZE, fmt, l);
  assert(res >= 0 && res < BUF_SIZE);
  va_end(l);
  return buffer;
}
