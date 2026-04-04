// Copyright 2016 Eliot Courtney.
#ifndef MODEL_BASE_H_
#define MODEL_BASE_H_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>

#include "util/error.h"

namespace mrna {

using Base = uint8_t;
using BaseMask = uint32_t;

// Used to index into RNA primary and secondary structures.
// Index stores both sequence positions and sequence lengths, so INVALID_INDEX
// is reserved as the sentinel value and the largest valid size/index is
// INVALID_INDEX - 1.
constexpr int INDEX_BITS = MRNA_INDEX_BITS;

constexpr uint8_t INVALID_UINT8 = std::numeric_limits<uint8_t>::max();
constexpr uint16_t INVALID_UINT16 = std::numeric_limits<uint16_t>::max();
constexpr uint32_t INVALID_UINT32 = std::numeric_limits<uint32_t>::max();

#if MRNA_INDEX_BITS == 8
using Index = uint8_t;
constexpr Index INVALID_INDEX = static_cast<Index>(INVALID_UINT8);
#elif MRNA_INDEX_BITS == 16
using Index = uint16_t;
constexpr Index INVALID_INDEX = static_cast<Index>(INVALID_UINT16);
#elif MRNA_INDEX_BITS == 32
using Index = int32_t;
constexpr Index INVALID_INDEX = std::numeric_limits<int32_t>::max();
#else
#error "Unsupported MRNA_INDEX_BITS"
#endif

constexpr std::size_t MAX_RNA_SIZE = std::size_t(INVALID_INDEX) - 1;

constexpr void VerifyRnaSize(std::size_t size) {
  verify(size <= MAX_RNA_SIZE, "RNA too long - recompile with larger INDEX_BITS");
}

// Don't ever change these values.
constexpr Base A = 0, C = 1, G = 2, U = 3, MAX_BASE = 4, MIN_BASE = 0;

constexpr BaseMask A_b = 1 << A, C_b = 1 << C, G_b = 1 << G, U_b = 1 << U;

constexpr bool IsBase(Base b) { return b <= 3; }

constexpr bool IsPairOf(Base a, Base b, BaseMask basesA, BaseMask basesB) {
  return (basesA & (1 << a)) && (basesB & (1 << b));
}

constexpr bool IsUnorderedOf(Base a, Base b, BaseMask basesA, BaseMask basesB) {
  return IsPairOf(a, b, basesA, basesB) || IsPairOf(b, a, basesA, basesB);
}

constexpr bool IsPair(Base a, Base b) {
  const BaseMask combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (G_b | C_b) || combined == (A_b | U_b);
}

// Returns if the opening and closing base pairs are considered to be continous.
// This happens if they form a stack or a size one bulge loop.
constexpr bool IsContinuous(int ost, int oen, int ist, int ien) {
  assert(ost < oen && ist < ien && ost < ist && oen > ien && ost >= 0 && oen >= 0);
  return (ost + 1 == ist && oen - 1 == ien) || (ost + 1 == ist && oen - 2 == ien) ||
      (ost + 2 == ist && oen - 1 == ien);
}

constexpr bool IsAuPair(Base a, Base b) { return (a == A && b == U) || (a == U && b == A); }

constexpr bool IsGuPair(Base a, Base b) { return (a == G && b == U) || (a == U && b == G); }

constexpr bool IsWcPair(Base a, Base b) {
  const BaseMask combined = (1 << a) | (1 << b);
  return combined == (A_b | U_b) || combined == (G_b | C_b);
}

constexpr bool IsGu(Base b) { return b == G || b == U; }

// Returns the Watson-Crick complement of the given base.
constexpr Base WcPair(Base b) { return b ^ 3; }

// Returns the GU complement of the given base.
constexpr Base GuPair(Base b) {
  assert(IsGu(b));
  return b ^ 1;
}

std::optional<Base> CharToBase(char c);

std::optional<char> BaseToChar(Base b);

}  // namespace mrna

#endif  // MODEL_BASE_H_
