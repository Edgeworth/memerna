// Copyright 2026 Eliot Courtney.
#ifndef UTIL_LOG_H_
#define UTIL_LOG_H_

#include <fmt/core.h>

#include <cstdio>
#include <utility>

#ifdef MRNA_ENABLE_LOGGING
#include <spdlog/spdlog.h>
#endif

namespace mrna {

enum class LogLevel {
  TRACE,
  DEBUG,
  INFO,
  WARN,
  ERR,
  CRITICAL,
  OFF,
};

void InitLog();
void SetLogLevel(LogLevel level);
[[nodiscard]] LogLevel GetLogLevel();

template <typename... Args>
inline void IgnoreUnused(const Args&...) {}

template <typename... Args>
inline void logdebug(fmt::format_string<Args...> fmt_str, Args&&... args) {
#ifdef MRNA_ENABLE_LOGGING
  spdlog::debug(fmt_str, std::forward<Args>(args)...);
#else
  (void)fmt_str;
  IgnoreUnused(args...);
#endif
}

template <typename... Args>
inline void loginfo(fmt::format_string<Args...> fmt_str, Args&&... args) {
#ifdef MRNA_ENABLE_LOGGING
  spdlog::info(fmt_str, std::forward<Args>(args)...);
#else
  (void)fmt_str;
  IgnoreUnused(args...);
#endif
}

template <typename... Args>
inline void logwarn(fmt::format_string<Args...> fmt_str, Args&&... args) {
#ifdef MRNA_ENABLE_LOGGING
  spdlog::warn(fmt_str, std::forward<Args>(args)...);
#else
  (void)fmt_str;
  IgnoreUnused(args...);
#endif
}

template <typename... Args>
inline void logcritical(fmt::format_string<Args...> fmt_str, Args&&... args) {
#ifdef MRNA_ENABLE_LOGGING
  spdlog::critical(fmt_str, std::forward<Args>(args)...);
#else
  fmt::print(stderr, fmt_str, std::forward<Args>(args)...);
#endif
}

}  // namespace mrna

#endif  // UTIL_LOG_H_
