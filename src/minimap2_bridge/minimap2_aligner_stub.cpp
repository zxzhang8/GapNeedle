#include "gapneedle/aligner.hpp"

#include <filesystem>
#include <stdexcept>

namespace gapneedle {

AlignmentResult Minimap2Aligner::align(const AlignmentRequest& request) const {
  if (request.reuseExisting && !request.outputPafPath.empty() &&
      std::filesystem::exists(request.outputPafPath)) {
    return AlignmentResult{request.outputPafPath, true, {"reused existing paf"}};
  }
  throw std::runtime_error(
      "minimap2 source is not integrated yet. Add minimap2 sources under third_party/minimap2 "
      "or build with GAPNEEDLE_USE_MINIMAP2=OFF and provide existing PAF for stitching.");
}

}  // namespace gapneedle
