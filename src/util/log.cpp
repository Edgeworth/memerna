// Copyright 2026 Eliot Courtney.
#include "util/log.h"

#include <cstdlib>
#include <utility>

#ifdef MRNA_ENABLE_LOGGING
#include "spdlog/cfg/env.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "util/version.h"
#endif

namespace mrna {

namespace {

#ifdef MRNA_ENABLE_LOGGING
spdlog::level::level_enum ToSpdlogLevel(LogLevel level) {
  switch (level) {
  case LogLevel::TRACE: return spdlog::level::trace;
  case LogLevel::DEBUG: return spdlog::level::debug;
  case LogLevel::INFO: return spdlog::level::info;
  case LogLevel::WARN: return spdlog::level::warn;
  case LogLevel::ERR: return spdlog::level::err;
  case LogLevel::CRITICAL: return spdlog::level::critical;
  case LogLevel::OFF: return spdlog::level::off;
  }
  std::abort();
}

LogLevel FromSpdlogLevel(spdlog::level::level_enum level) {
  switch (level) {
  case spdlog::level::trace: return LogLevel::TRACE;
  case spdlog::level::debug: return LogLevel::DEBUG;
  case spdlog::level::info: return LogLevel::INFO;
  case spdlog::level::warn: return LogLevel::WARN;
  case spdlog::level::err: return LogLevel::ERR;
  case spdlog::level::critical: return LogLevel::CRITICAL;
  case spdlog::level::off: return LogLevel::OFF;
  case spdlog::level::n_levels: break;
  }
  std::abort();
}
#else
LogLevel level = LogLevel::INFO;
#endif

}  // namespace

void InitLog() {
#ifdef MRNA_ENABLE_LOGGING
  auto logger = spdlog::get("stderr");
  if (!logger) logger = spdlog::stderr_color_mt("stderr");
  spdlog::set_default_logger(std::move(logger));
  spdlog::cfg::load_env_levels();
  loginfo("memerna version {}", VERSION.ToString());
#endif
}

void SetLogLevel(LogLevel new_level) {
#ifdef MRNA_ENABLE_LOGGING
  spdlog::set_level(ToSpdlogLevel(new_level));
#else
  level = new_level;
#endif
}

LogLevel GetLogLevel() {
#ifdef MRNA_ENABLE_LOGGING
  return FromSpdlogLevel(spdlog::get_level());
#else
  return level;
#endif
}

}  // namespace mrna
