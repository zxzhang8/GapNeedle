#pragma once

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace gapneedle {

struct AlignmentRequest {
  std::string targetFasta;
  std::string queryFasta;
  std::string targetSeq;
  std::string querySeq;
  std::string preset{"asm10"};
  int threads{4};
  bool reverseTarget{false};
  bool reverseQuery{false};
  bool reuseExisting{true};
  std::string outputPafPath;
};

struct AlignmentRecord {
  std::string qName;
  int qLen{0};
  int qStart{0};
  int qEnd{0};
  char strand{'+'};
  std::string tName;
  int tLen{0};
  int tStart{0};
  int tEnd{0};
  int matches{0};
  int alnLen{0};
  int mapq{0};
  std::vector<std::string> extras;
};

struct AlignmentResult {
  std::string pafPath;
  bool skipped{false};
  std::vector<std::string> warnings;
};

struct Segment {
  std::string source;   // t, q, x1...
  std::string seqName;
  int start{0};
  int end{0};
  bool reverse{false};
};

struct BreakpointSummary {
  int index{0};
  bool leftFlankMatch{false};
  bool rightFlankMatch{false};
  std::string preview;
};

struct StitchRequest {
  std::string targetFasta;
  std::string queryFasta;
  std::unordered_map<std::string, std::string> extraFastaBySource;
  std::vector<Segment> segments;
  int contextBp{200};
  std::string outputFastaPath;
  std::string outputSeqName;
};

struct StitchResult {
  std::string outputFastaPath;
  std::string outputLogPath;
  std::size_t mergedLength{0};
  std::vector<BreakpointSummary> breakpoints;
};

struct MappingResult {
  std::optional<int> tPos;
  std::string reason;
  int qPos{0};
  std::optional<int> qPosOriented;
  char op{0};
  int opLen{0};
  int opOffset{0};
  std::unordered_map<char, int> countsBefore;
  std::unordered_map<char, int> countsTotal;
  int qConsumedBefore{0};
  int tConsumedBefore{0};
};

}  // namespace gapneedle
