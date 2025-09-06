// Copyright 2023 Eliot Courtney.
#include "util/error.h"

#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <cstdlib>
#include <exception>
#include <ios>

#include "spdlog/sinks/stdout_color_sinks.h"
#include "util/version.h"

#ifdef USE_BOOST
#include <boost/stacktrace.hpp>
#endif

namespace {
void terminate_handler() {
  std::exception_ptr exptr = std::current_exception();
  if (exptr) {
    try {
      std::rethrow_exception(exptr);
    } catch (std::exception &ex) {
      fmt::print(stderr, "terminated due to exception: {}\n", ex.what());
    } catch (...) {
      fmt::print(stderr, "terminated due to unknown exception\n");
    }
  } else {
    fmt::print(stderr, "terminated due to unknown reason\n");
  }

#ifdef USE_BOOST
  auto stack = boost::stacktrace::to_string(boost::stacktrace::stacktrace());
  fmt::print(stderr, "stack trace:\n{}\n", stack);
#endif

  std::abort();
}

}  // namespace

namespace mrna {

void InitProgram() {
  std::ios_base::sync_with_stdio(false);
  std::set_terminate(terminate_handler);
  auto logger = spdlog::stderr_color_mt("stderr");
  spdlog::set_default_logger(logger);
  spdlog::info("memerna version {}", VERSION.ToString());
}

}  // namespace mrna
