// Copyright 2021 Eliot Courtney.
#ifndef UTIL_FLOAT_H_
#define UTIL_FLOAT_H_

#include <algorithm>
#include <cmath>

#ifdef USE_MPFR
#include <fmt/ostream.h>

#include <boost/multiprecision/mpfr.hpp>

#include "util/string.h"

namespace boost {

inline void throw_exception(const std::exception& e) {
  fmt::print(stderr, "Boost exception: {}\n", e.what());
  std::abort();
}

}  // namespace boost
#endif  // USE_MPFR

namespace mrna {

constexpr int powi(int base, int exp) {
  int result = 1;
  for (int i = 0; i < exp; ++i) result *= base;
  return result;
}

#ifdef USE_MPFR

using flt =
    boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<FLOAT_PRECISION>>;

#define _STRHELP(x) STRINGIFY(x)
#define FLTFMT "{:." _STRHELP(FLOAT_PRECISION) "f}"

// Stringify to ensure
#define FLT(x) static_cast<flt>(#x)

#else

#if FLOAT_PRECISION == 6
using flt = float;
#define FLTFMT "{:.6f}"
#elif FLOAT_PRECISION == 15
using flt = double;
#define FLTFMT "{:.15f}"
#elif FLOAT_PRECISION == 18
using flt = long double;
#define FLTFMT "{:.18f}"
#else
static_assert(false, "unknown float precision");
#endif

// Redefinitions to use std:: for long double support, but also to keep
// MPFR happy.
inline flt fabs(flt v) { return std::abs(v); }
inline flt exp(flt v) { return std::exp(v); }

// Just pass through for built in types. Make literals long doubles.
#define FLT(x) static_cast<flt>(x##l)

#endif  // USE_MPFR

#if FLOAT_PRECISION == 6
inline const flt EP{1e-3};
#elif FLOAT_PRECISION == 15
inline const flt EP{1e-7};
#elif FLOAT_PRECISION == 18
inline const flt EP{1e-9};
#else
inline const flt EP{1e-30};
#endif

inline bool abs_eq(flt a, flt b, flt ep = EP) { return fabs(a - b) < ep; }

inline bool rel_eq(flt a, flt b, flt rel = EP) {
  return fabs(a - b) <= rel * std::max(fabs(a), fabs(b));
}

inline bool absrel_eq(flt a, flt b, flt ep = EP) { return abs_eq(a, b, ep) || rel_eq(a, b, ep); }

}  // namespace mrna

#ifdef USE_MPFR

#include <fmt/format.h>

template <typename T>
concept FltConvertible = boost::multiprecision::is_number<T>::value ||
    boost::multiprecision::is_number_expression<T>::value;

namespace fmt {
template <FltConvertible T>
struct formatter<T> {
 private:
  detail::dynamic_format_specs<char> specs_;

 public:
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    auto end = parse_format_specs(ctx.begin(), ctx.end(), specs_, ctx, detail::type::float_type);
    detail::parse_float_type_spec(specs_);
    return end;
  }

  auto format(const T& c, format_context& ctx) const -> auto {
    std::stringstream out;
    int precision = specs_.precision;
    detail::handle_dynamic_spec<detail::precision_checker>(precision, specs_.precision_ref, ctx);
    if (precision > 0) out << std::setprecision(precision);
    if (specs_.width > 0) out << std::setw(specs_.width);
    if (specs_.type == presentation_type::fixed) out << std::fixed;
    if (specs_.type == presentation_type::exp) out << std::scientific;
    out << c;
    return format_to(ctx.out(), "{}", out.str());
  }
};
}  // namespace fmt

#endif  // USE_MPFR

#endif  // UTIL_FLOAT_H_
