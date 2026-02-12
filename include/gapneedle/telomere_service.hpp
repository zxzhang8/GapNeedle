#pragma once

#include <string>
#include <utility>

namespace gapneedle {

std::pair<bool, bool> checkTelomere(const std::string& fastaPath,
                                    const std::string& seqName,
                                    int window = 1000000,
                                    const std::string& motif = "CCCTAA",
                                    int minRepeats = 15);

}  // namespace gapneedle
