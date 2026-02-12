#include "gapneedle/facade.hpp"

#include "gapneedle/fasta_io.hpp"

namespace gapneedle {

GapNeedleFacade::GapNeedleFacade() = default;

AlignmentResult GapNeedleFacade::align(const AlignmentRequest& request) const {
  return aligner_.align(request);
}

StitchResult GapNeedleFacade::stitch(const StitchRequest& request) const {
  return stitchService_.stitch(request);
}

std::vector<std::tuple<std::string, int, int>> GapNeedleFacade::scanGaps(const std::string& fastaPath,
                                                                          int minGap) const {
  std::vector<std::tuple<std::string, int, int>> gaps;
  auto records = readFasta(fastaPath);
  for (const auto& [name, seq] : records) {
    int start = -1;
    for (int i = 0; i < static_cast<int>(seq.size()); ++i) {
      if (seq[i] == 'N') {
        if (start < 0) {
          start = i;
        }
      } else {
        if (start >= 0 && i - start >= minGap) {
          gaps.emplace_back(name, start, i);
        }
        start = -1;
      }
    }
    if (start >= 0 && static_cast<int>(seq.size()) - start >= minGap) {
      gaps.emplace_back(name, start, static_cast<int>(seq.size()));
    }
  }
  return gaps;
}

}  // namespace gapneedle
