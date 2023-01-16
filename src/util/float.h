// Copyright 2021 Eliot Courtney.
#ifndef UTIL_FLOAT_H_
#define UTIL_FLOAT_H_

#include <algorithm>
#include <cmath>
#include <cstdint>

#ifdef USE_MPFR
#include <fmt/ostream.h>

#include <boost/multiprecision/mpfr.hpp>

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
using flt = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<FLOAT_BITS>>;
#else
#if FLOAT_BITS == 32
using flt = float;
#elif FLOAT_BITS == 64
using flt = double;
#elif FLOAT_BITS == 80
using flt = long double;
#else
static_assert(false, "unknown float bits")
#endif

// Redefinitions to use std:: for long double support, but also to keep
// MPFR happy.
inline flt fabs(flt v) { return std::abs(v); }
inline flt exp(flt v) { return std::exp(v); }

#endif  // USE_MPFR

#if FLOAT_BITS == 32
inline const flt EP{1e-3};
#elif FLOAT_BITS == 64
inline const flt EP{1e-6};
#elif FLOAT_BITS == 80
inline const flt EP{1e-8};
#else
inline const flt EP{1e-30};
#endif

inline bool abs_eq(flt a, flt b, flt ep = EP) { return fabs(a - b) < ep; }

inline bool rel_eq(flt a, flt b, flt rel = EP) {
  return fabs(a - b) <= rel * std::max(fabs(a), fabs(b));
}

}  // namespace mrna

#ifdef USE_MPFR

namespace fmt {
template <>
struct formatter<::mrna::flt> {
 private:
  detail::dynamic_format_specs<char> specs_;

 public:
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    using handler_type = detail::dynamic_specs_handler<format_parse_context>;
    auto type = detail::type_constant<::mrna::flt, char>::value;
    auto checker = detail::specs_checker<handler_type>(handler_type(specs_, ctx), type);
    auto end = detail::parse_format_specs(ctx.begin(), ctx.end(), checker);
    return end;
  }

  template <class FormatContext>
  auto format(const ::mrna::flt& c, FormatContext& ctx) -> decltype(ctx.out()) {
    std::stringstream out;
    detail::handle_dynamic_spec<detail::precision_checker>(
        specs_.precision, specs_.precision_ref, ctx);
    if (specs_.precision > 0) out << std::setprecision(specs_.precision);
    if (specs_.width > 0) out << std::setw(specs_.width);
    if (specs_.type == presentation_type::fixed_lower) out << std::fixed;
    if (specs_.type == presentation_type::exp_lower) out << std::scientific;
    out << c;
    return format_to(ctx.out(), "{}", out.str());
  }
};
}  // namespace fmt

#endif  // USE_MPFR

#endif  // UTIL_FLOAT_H_
