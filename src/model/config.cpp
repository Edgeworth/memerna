// Copyright 2021 E.
// Copyright 2021 E.
#include "model/config.h"

namespace mrna {

context_opt_t ContextOptionsFromArgParse(const ArgParse& args) {
  context_opt_t options;
  const auto dp_alg = args.GetOption("dp-alg");
  if (dp_alg == "0") {
    options.table_alg = context_opt_t::TableAlg::ZERO;
  } else if (dp_alg == "1") {
    options.table_alg = context_opt_t::TableAlg::ONE;
  } else if (dp_alg == "2") {
    options.table_alg = context_opt_t::TableAlg::TWO;
  } else if (dp_alg == "3") {
    options.table_alg = context_opt_t::TableAlg::THREE;
  } else if (dp_alg == "brute") {
    options.table_alg = context_opt_t::TableAlg::BRUTE;
  } else {
    error("unknown fold option");
  }
  const auto subopt_alg = args.GetOption("subopt-alg");
  if (subopt_alg == "0") {
    options.suboptimal_alg = context_opt_t::SuboptimalAlg::ZERO;
  } else if (subopt_alg == "1") {
    options.suboptimal_alg = context_opt_t::SuboptimalAlg::ONE;
  } else if (subopt_alg == "brute") {
    options.suboptimal_alg = context_opt_t::SuboptimalAlg::BRUTE;
  } else {
    error("unknown suboptimal option");
  }
  const auto part_alg = args.GetOption("part-alg");
  if (part_alg == "0") {
    options.partition_alg = context_opt_t::PartitionAlg::ZERO;
  } else if (part_alg == "1") {
    options.partition_alg = context_opt_t::PartitionAlg::ONE;
  } else if (part_alg == "brute") {
    options.partition_alg = context_opt_t::PartitionAlg::BRUTE;
  } else {
    error("unknown partition option");
  }
  return options;
}

}  // namespace mrna
