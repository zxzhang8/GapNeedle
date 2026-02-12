#pragma once

#include "gapneedle/types.hpp"

#include <string>
#include <vector>

namespace gapneedle {

std::vector<AlignmentRecord> parsePaf(const std::string& path,
                                      const std::string& targetSeq,
                                      const std::string& querySeq);
std::vector<AlignmentRecord> suggestOverlaps(const std::string& path,
                                             const std::string& targetSeq,
                                             const std::string& querySeq,
                                             int limit = 10);

}  // namespace gapneedle
