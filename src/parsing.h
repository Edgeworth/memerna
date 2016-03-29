#ifndef MEMERNA_PARSING_H
#define MEMERNA_PARSING_H

#include <string>
#include <stack>
#include "base.h"

namespace memerna {
namespace parsing {

rna_t ParseRnaFromString(const std::string& s);

folded_rna_t ParseViennaRna(const std::string& rna_str, const std::string& pairs_str);

void ParseStackingEnergiesFromFile(const std::string& filename);

}
}

#endif //MEMERNA_PARSING_H
