// Copyright 2021 Eliot Courtney.
#ifndef UTIL_FLOAT_H_
#define UTIL_FLOAT_H_

#include <algorithm>
#include <cmath>
#include <string>

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

inline flt ParseFlt(const std::string& str) { return flt(str); }

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
inline flt round(flt v) { return std::round(v); }
inline flt isnan(flt v) { return std::isnan(v); }

// Just pass through for built in types. Make literals long doubles.
#define FLT(x) static_cast<flt>(x##l)

inline flt ParseFlt(const std::string& str) { return static_cast<flt>(std::stold(str)); }

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
    return parse_format_specs(ctx.begin(), ctx.end(), specs_, ctx, detail::type::float_type);
  }

  auto format(const T& val, format_context& ctx) const -> auto {
    auto specs = specs_;

    detail::handle_dynamic_spec(specs.dynamic_width(), specs.width, specs.width_ref, ctx);
    detail::handle_dynamic_spec(
        specs.dynamic_precision(), specs.precision, specs.precision_ref, ctx);

    if (!boost::multiprecision::isfinite(val)) {
      sign s = boost::multiprecision::signbit(val) ? sign::minus : specs.sign();
      return detail::write_nonfinite<char>(ctx.out(), boost::multiprecision::isnan(val), specs, s);
    }

    int precision = specs.precision >= 0 ? specs.precision : 6;
    std::ios_base::fmtflags flags = std::ios_base::dec;
    if (specs.type() == presentation_type::fixed)
      flags |= std::ios_base::fixed;
    else if (specs.type() == presentation_type::exp)
      flags |= std::ios_base::scientific;

    std::string num = mrna::flt(val).str(precision, flags);

    const char* prefix = "";
    if (boost::multiprecision::signbit(val)) {
      prefix = "-";
      if (!num.empty() && num[0] == '-') num.erase(0, 1);
    } else {
      switch (specs.sign()) {
      case sign::plus: prefix = "+"; break;
      case sign::space: prefix = " "; break;
      default: break;
      }
    }

    size_t num_size = num.size();
    size_t prefix_size = std::strlen(prefix);
    size_t total_size = num_size + prefix_size;

    using iterator = decltype(ctx.out());
    return detail::write_padded<char>(ctx.out(), specs, total_size, total_size, [&](iterator it) {
      if (prefix_size != 0) it = std::copy_n(prefix, prefix_size, it);
      it = std::copy_n(num.data(), num_size, it);
      return it;
    });
  }
};

}  // namespace fmt

#endif  // USE_MPFR

#endif  // UTIL_FLOAT_H_
