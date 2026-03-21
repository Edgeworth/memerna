// Copyright 2025 Eliot Courtney.
#ifndef UTIL_VERSION_H_
#define UTIL_VERSION_H_

#include <fmt/core.h>

#include <string>

namespace mrna {

struct VersionInfo {
  const unsigned int major;
  const unsigned int minor;
  const unsigned int patch;
  const std::string pre_release;
  const std::string build_metadata;

  [[nodiscard]] std::string ToString() const {
    return fmt::format("{}.{}.{}{}{}", major, minor, patch,
        pre_release.empty() ? "" : "-" + std::string(pre_release),
        build_metadata.empty() ? "" : "+" + std::string(build_metadata));
  }
};

inline constexpr VersionInfo VERSION = {.major = MRNA_VERSION_MAJOR,
    .minor = MRNA_VERSION_MINOR,
    .patch = MRNA_VERSION_PATCH,
    .pre_release = "",
    .build_metadata = ""};

}  // namespace mrna

#endif  // UTIL_VERSION_H_
