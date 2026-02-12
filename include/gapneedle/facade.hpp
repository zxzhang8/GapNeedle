#pragma once

#include "gapneedle/aligner.hpp"
#include "gapneedle/stitch_service.hpp"
#include "gapneedle/types.hpp"

#include <tuple>
#include <vector>

namespace gapneedle {

class GapNeedleFacade {
 public:
  GapNeedleFacade();

  AlignmentResult align(const AlignmentRequest& request) const;
  StitchResult stitch(const StitchRequest& request) const;
  std::vector<std::tuple<std::string, int, int>> scanGaps(const std::string& fastaPath,
                                                          int minGap = 10) const;

 private:
  Minimap2Aligner aligner_;
  StitchService stitchService_;
};

}  // namespace gapneedle
