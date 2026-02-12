#include "gapneedle/telomere_service.hpp"

#include "gapneedle/fasta_io.hpp"

#include <stdexcept>

namespace gapneedle {

namespace {

bool hasConsecutiveMotif(const std::string& seq, const std::string& motif, int minRepeats) {
  if (motif.empty() || minRepeats <= 0) {
    return false;
  }
  int count = 0;
  std::size_t pos = 0;
  while (pos + motif.size() <= seq.size()) {
    if (seq.compare(pos, motif.size(), motif) == 0) {
      ++count;
      if (count >= minRepeats) {
        return true;
      }
      pos += motif.size();
    } else {
      count = 0;
      ++pos;
    }
  }
  return false;
}

}  // namespace

std::pair<bool, bool> checkTelomere(const std::string& fastaPath,
                                    const std::string& seqName,
                                    int window,
                                    const std::string& motif,
                                    int minRepeats) {
  const auto seqs = readFastaSelected(fastaPath, {seqName});
  auto it = seqs.find(seqName);
  if (it == seqs.end()) {
    throw std::runtime_error("Sequence not found: " + seqName);
  }
  const std::string& seq = it->second;
  const int w = std::min<int>(window, static_cast<int>(seq.size()));
  const std::string left = seq.substr(0, w);
  const std::string right = seq.substr(seq.size() - w);

  const std::string motifUpper = motif;
  const std::string rc = reverseComplement(motifUpper);

  const bool leftHas = hasConsecutiveMotif(left, motifUpper, minRepeats) ||
                       hasConsecutiveMotif(left, rc, minRepeats);
  const bool rightHas = hasConsecutiveMotif(right, motifUpper, minRepeats) ||
                        hasConsecutiveMotif(right, rc, minRepeats);
  return {leftHas, rightHas};
}

}  // namespace gapneedle
