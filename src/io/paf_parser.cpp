#include "gapneedle/paf.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace gapneedle {

std::vector<AlignmentRecord> parsePaf(const std::string& path,
                                      const std::string& targetSeq,
                                      const std::string& querySeq) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open PAF: " + path);
  }

  std::vector<AlignmentRecord> out;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, '\t')) {
      parts.push_back(item);
    }
    if (parts.size() < 12) {
      continue;
    }
    if (parts[0] != querySeq || parts[5] != targetSeq) {
      continue;
    }

    AlignmentRecord r;
    r.qName = parts[0];
    r.qLen = std::stoi(parts[1]);
    r.qStart = std::stoi(parts[2]);
    r.qEnd = std::stoi(parts[3]);
    r.strand = parts[4].empty() ? '+' : parts[4][0];
    r.tName = parts[5];
    r.tLen = std::stoi(parts[6]);
    r.tStart = std::stoi(parts[7]);
    r.tEnd = std::stoi(parts[8]);
    r.matches = std::stoi(parts[9]);
    r.alnLen = std::stoi(parts[10]);
    r.mapq = std::stoi(parts[11]);
    for (std::size_t i = 12; i < parts.size(); ++i) {
      r.extras.push_back(parts[i]);
    }
    out.push_back(std::move(r));
  }
  return out;
}

std::vector<AlignmentRecord> suggestOverlaps(const std::string& path,
                                             const std::string& targetSeq,
                                             const std::string& querySeq,
                                             int limit) {
  auto records = parsePaf(path, targetSeq, querySeq);
  std::sort(records.begin(), records.end(), [](const AlignmentRecord& a, const AlignmentRecord& b) {
    int oa = std::min(a.qEnd - a.qStart, a.tEnd - a.tStart);
    int ob = std::min(b.qEnd - b.qStart, b.tEnd - b.tStart);
    return oa > ob;
  });
  if (limit > 0 && static_cast<int>(records.size()) > limit) {
    records.resize(static_cast<std::size_t>(limit));
  }
  return records;
}

}  // namespace gapneedle
