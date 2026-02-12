#pragma once

#include "gapneedle/types.hpp"

namespace gapneedle {

class IAligner {
 public:
  virtual ~IAligner() = default;
  virtual AlignmentResult align(const AlignmentRequest& request) const = 0;
};

class Minimap2Aligner final : public IAligner {
 public:
  AlignmentResult align(const AlignmentRequest& request) const override;
};

}  // namespace gapneedle
