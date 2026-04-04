// Copyright 2022 Eliot Courtney.
#include "tests/init.h"

#include <string>
#include <tuple>

#include "api/ctx/backend.h"
#include "api/energy/energy_cfg.h"

namespace mrna {

std::tuple<Primary, Secondary> kNNDBHairpin1 = ParseSeqDb("CACAAAAAAAUGUG", "((((......))))");
std::tuple<Primary, Secondary> kNNDBHairpin2 = ParseSeqDb("CACAGGAAGUGUG", "((((.....))))");
std::tuple<Primary, Secondary> kNNDBHairpin3 = ParseSeqDb("CACCCGAGGGUG", "((((....))))");
std::tuple<Primary, Secondary> kNNDBHairpin4 = ParseSeqDb("CACACCCCCCUGUG", "((((......))))");
std::tuple<Primary, Secondary> kNNDBHairpin5 = ParseSeqDb("CGGGGGAAGUCCG", "((((.....))))");
std::tuple<Primary, Secondary> kNNDBBulge1 = ParseSeqDb("GCCCGAAACGGC", "(((.(...))))");
std::tuple<Primary, Secondary> kNNDBBulge2 = ParseSeqDb("GAACAGAAACUC", "((...(...)))");
std::tuple<Primary, Secondary> kNNDBInternal2x3 =
    ParseSeqDb("CAGACGAAACGGAGUG", "((..((...))...))");
std::tuple<Primary, Secondary> kNNDBInternal1x5 =
    ParseSeqDb("CAGCGAAACGGAAAGUG", "((.((...)).....))");
std::tuple<Primary, Secondary> kNNDBInternal2x2 = ParseSeqDb("CAGACGAAACGGAUG", "((..((...))..))");
std::tuple<Primary, Secondary> kBulge1 = ParseSeqDb("GCUCGAAACAGC", "(((.(...))))");
std::tuple<Primary, Secondary> kInternal1 = ParseSeqDb("AGAGAAACAAAU", "(..(...)...)");
std::optional<std::tuple<Primary, Secondary>> k16sHSapiens3;

void InitTest(const std::string& data_dir) {
  // 315 nt sequence - only available when MAX_RNA_SIZE is large enough.
  constexpr const char* kK16sHSapiens3Seq =
      "AAGGACCUGGCGGUGCUUCAUAUCCCUCUAGAGGAGCCUGUUCUGUAAUCGAUAAACCCCGAUCAACCUCACCACCUCUUGCU"
      "CAGCCUAUAUACCGCCAUCUUCAGCAAACCCUGAUGAAGGCUACAAAGUAAGCGCAAGUACCCACGUAAAGACGUUAGGUCAA"
      "GGUGUAGCCCAUGAGGUGGCAAGAAAUGGGCUACAUUUUCUACCCCAGAAAACUACGAUAGCCCUUAUGAAACUUAAGGGUCG"
      "AAGGUGGAUUUAGCAGUAAACUAAGAGUAGAGUGCUUAGUUGAACAGGGCCCUGAAGCGCGUACAC";
  constexpr const char* kK16sHSapiens3Db =
      ".......(((((.(((((((...((..((((((.((((((((((...((((........))))........(((((((.......((.(("
      "((..((((((.(.(..((..(((((.....))).......))..)).....(((....)))...).).).)))...))))))))....))"
      ")))))..)).)))))))).)...((((.....)))).....(..(.(((((((.......))))))).)..).....))))).....((("
      "((((.........)))))))......))...)))))))))).)).";
  if (std::char_traits<char>::length(kK16sHSapiens3Seq) <= MAX_RNA_SIZE)
    k16sHSapiens3 = ParseSeqDb(kK16sHSapiens3Seq, kK16sHSapiens3Db);

  // All backends support T04 and T12.
  for (auto backend : EnumValues<BackendKind>()) {
    t04_ms.push_back(BackendFromBackendCfg(backend,
        BackendCfg{.energy_model = erg::EnergyModelKind::T04,
            .precision = MRNA_ENERGY_PRECISION,
            .data_src = data_dir}));
    t12_ms.push_back(BackendFromBackendCfg(backend,
        BackendCfg{.energy_model = erg::EnergyModelKind::T12,
            .precision = MRNA_ENERGY_PRECISION,
            .data_src = data_dir}));
  }

  // Only STACK supports T22.
  t22_ms.push_back(BackendFromBackendCfg(BackendKind::STACK,
      BackendCfg{.energy_model = erg::EnergyModelKind::T22,
          .precision = MRNA_ENERGY_PRECISION,
          .data_src = data_dir}));

  verify(t04_ms.size() == NUM_T04_MODELS, "t04_ms.size() == NUM_T04_MODELS");
  verify(t12_ms.size() == NUM_T12_MODELS, "t12_ms.size() == NUM_T12_MODELS");
  verify(t22_ms.size() == NUM_T22_MODELS, "t22_ms.size() == NUM_T22_MODELS");

  base_t04 = md::base::Model::FromBackendCfg(BackendKind::BASE,
      BackendCfg{.energy_model = erg::EnergyModelKind::T04,
          .precision = MRNA_ENERGY_PRECISION,
          .data_src = data_dir});
  base_ms.push_back(base_t04);
  while (base_ms.size() < NUM_TEST_MODELS)
    base_ms.push_back(md::base::Model::Random(kT04Cfg, base_ms.size()));

  baseopt_t04 = md::base::opt::Model::FromBackendCfg(BackendKind::BASEOPT,
      BackendCfg{.energy_model = erg::EnergyModelKind::T04,
          .precision = MRNA_ENERGY_PRECISION,
          .data_src = data_dir});
  baseopt_ms.push_back(baseopt_t04);
  while (baseopt_ms.size() < NUM_TEST_MODELS)
    baseopt_ms.push_back(md::base::opt::Model::Random(kT04Cfg, baseopt_ms.size()));

  stack_t04 = md::stack::Model::FromBackendCfg(BackendKind::STACK,
      BackendCfg{.energy_model = erg::EnergyModelKind::T04,
          .precision = MRNA_ENERGY_PRECISION,
          .data_src = data_dir});
  stack_ms.push_back(stack_t04);
  while (stack_ms.size() < NUM_TEST_MODELS)
    stack_ms.push_back(md::stack::Model::Random(kT04Cfg, stack_ms.size()));
}

}  // namespace mrna
