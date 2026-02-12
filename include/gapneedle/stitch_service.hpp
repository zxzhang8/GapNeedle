#pragma once

#include "gapneedle/types.hpp"

namespace gapneedle {

class StitchService {
 public:
  StitchResult stitch(const StitchRequest& request) const;
};

}  // namespace gapneedle
