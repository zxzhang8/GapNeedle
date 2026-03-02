#pragma once

#include "gapneedle/types.hpp"

namespace gapneedle {

class GuidedStitchService {
 public:
  GuidedSeedResult seedCandidates(const GuidedSeedRequest& request) const;
  GuidedStepResult nextCandidates(const GuidedStepRequest& request) const;
};

}  // namespace gapneedle
