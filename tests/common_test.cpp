// Copyright 2022 Eliot Courtney.
#include "common_test.h"

#include <variant>

#include "api/energy/energy.h"

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
std::tuple<Primary, Secondary> k16sHSapiens3 =
    ParseSeqDb("AAGGACCUGGCGGUGCUUCAUAUCCCUCUAGAGGAGCCUGUUCUGUAAUCGAUAAACCCCGAUCAACCUCACCACCUCUUGCU"
               "CAGCCUAUAUACCGCCAUCUUCAGCAAACCCUGAUGAAGGCUACAAAGUAAGCGCAAGUACCCACGUAAAGACGUUAGGUCAA"
               "GGUGUAGCCCAUGAGGUGGCAAGAAAUGGGCUACAUUUUCUACCCCAGAAAACUACGAUAGCCCUUAUGAAACUUAAGGGUCG"
               "AAGGUGGAUUUAGCAGUAAACUAAGAGUAGAGUGCUUAGUUGAACAGGGCCCUGAAGCGCGUACAC",
        ".......(((((.(((((((...((..((((((.((((((((((...((((........))))........(((((((.......((.(("
        "((..((((((.(.(..((..(((((.....))).......))..)).....(((....)))...).).).)))...))))))))....))"
        ")))))..)).)))))))).)...((((.....)))).....(..(.(((((((.......))))))).)..).....))))).....((("
        "((((.........)))))))......))...)))))))))).)).");

void InitTest(const std::string& data_dir) {
  int last_test_ems_idx = 0;
  int last_test_t04_ems_idx = 0;
  int last_test_t22_ems_idx = 0;

  // NEWMODEL: Add here.
#if ENERGY_PRECISION == 1
  mrna::t04p1 = mrna::md::t04::Model::FromModelPath(mrna::erg::ModelPath(data_dir, "t04p1"));

  mrna::test_ems[last_test_ems_idx++] = mrna::t04p1;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t04p1;

#elif ENERGY_PRECISION == 2
  mrna::t04p2 = mrna::md::t04::Model::FromModelPath(mrna::erg::ModelPath(data_dir, "t04p2"));
  mrna::t12p2 = mrna::md::t04::Model::FromModelPath(mrna::erg::ModelPath(data_dir, "t12p2"));
  mrna::t22p2 = mrna::md::t22::Model::FromModelPath(mrna::erg::ModelPath(data_dir, "t22p2"));

  mrna::test_ems[last_test_ems_idx++] = mrna::t04p2;
  mrna::test_ems[last_test_ems_idx++] = mrna::t12p2;
  mrna::test_ems[last_test_ems_idx++] = mrna::t22p2;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t04p2;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t12p2;
  mrna::test_t22_ems[last_test_t22_ems_idx++] = mrna::t22p2;

#endif

  for (int i = last_test_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_ems[i] = mrna::erg::Random(mrna::erg::ModelKind::T04_LIKE, i);

  for (int i = last_test_t04_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_t04_ems[i] = mrna::md::t04::Model::Random(i);

  for (int i = last_test_t22_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_t22_ems[i] = mrna::md::t22::Model::Random(i);
}

}  // namespace mrna
