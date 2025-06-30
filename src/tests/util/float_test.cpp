// Copyright 2025 Eliot Courtney.
#include "util/float.h"

#include <fmt/format.h>

#include <limits>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace mrna {

TEST(FloatTest, RoundTripFiniteValues) {
  const std::vector<flt> test_values = {
      FLT(0.0),
      FLT(1.0),
      FLT(-1.0),
      FLT(12345.6789),
      FLT(-9876.54321),
      FLT(4.548948615631896198153063518916351),
      FLT(1.2345e-20),
      FLT(-5.4321e25),
  };

  for (const auto& val : test_values) {
    std::string str = fmt::format(FLTFMT, val);
    flt parsed = ParseFlt(str);
    EXPECT_TRUE(absrel_eq(val, parsed));
  }
}

TEST(FloatTest, RoundTripSpecialValues) {
  if (std::numeric_limits<flt>::has_infinity) {
    const flt pos_inf = std::numeric_limits<flt>::infinity();
    std::string inf_str = fmt::format("{}", pos_inf);
    flt parsed_inf = ParseFlt(inf_str);
    EXPECT_EQ(pos_inf, parsed_inf);

    const flt neg_inf = -std::numeric_limits<flt>::infinity();
    std::string neg_inf_str = fmt::format("{}", neg_inf);
    flt parsed_neg_inf = ParseFlt(neg_inf_str);
    EXPECT_EQ(neg_inf, parsed_neg_inf);
  }

  if (std::numeric_limits<flt>::has_quiet_NaN) {
    const flt nan_val = std::numeric_limits<flt>::quiet_NaN();
    std::string nan_str = fmt::format("{}", nan_val);
    flt parsed_nan = ParseFlt(nan_str);
    EXPECT_TRUE(isnan(parsed_nan));
  }
}

}  // namespace mrna
