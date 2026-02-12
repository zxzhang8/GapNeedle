#include "gapneedle/stitch_service.hpp"

#include "gapneedle/fasta_io.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace gapneedle {

namespace {

std::string previewJunction(const std::string& left, const std::string& right, int ctx) {
  const int c = std::max(0, ctx);
  const std::string l = left.substr(left.size() > static_cast<std::size_t>(c) ? left.size() - c : 0);
  const std::string r = right.substr(0, std::min<int>(c, static_cast<int>(right.size())));
  return l + "|" + r;
}

bool compareSuffix(const std::string& a, const std::string& b, int n) {
  if (a.empty() || b.empty() || n <= 0) {
    return false;
  }
  n = std::min(n, std::min(static_cast<int>(a.size()), static_cast<int>(b.size())));
  return a.substr(a.size() - n) == b.substr(b.size() - n);
}

bool comparePrefix(const std::string& a, const std::string& b, int n) {
  if (a.empty() || b.empty() || n <= 0) {
    return false;
  }
  n = std::min(n, std::min(static_cast<int>(a.size()), static_cast<int>(b.size())));
  return a.substr(0, n) == b.substr(0, n);
}

std::string toJsonString(const std::string& raw) {
  std::string out;
  out.reserve(raw.size() + 8);
  for (char ch : raw) {
    switch (ch) {
      case '\\': out += "\\\\"; break;
      case '"': out += "\\\""; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default: out += ch; break;
    }
  }
  return out;
}

}  // namespace

StitchResult StitchService::stitch(const StitchRequest& request) const {
  if (request.segments.empty()) {
    throw std::runtime_error("stitch request has no segments");
  }
  if (request.outputFastaPath.empty()) {
    throw std::runtime_error("outputFastaPath is required");
  }

  const auto target = readFasta(request.targetFasta);
  const auto query = readFasta(request.queryFasta);

  std::unordered_map<std::string, FastaMap> extras;
  for (const auto& [src, path] : request.extraFastaBySource) {
    extras[src] = readFasta(path);
  }

  std::vector<std::string> pieceSeqs;
  pieceSeqs.reserve(request.segments.size());
  for (const auto& seg : request.segments) {
    const FastaMap* sourceMap = nullptr;
    if (seg.source == "t") {
      sourceMap = &target;
    } else if (seg.source == "q") {
      sourceMap = &query;
    } else {
      auto it = extras.find(seg.source);
      if (it == extras.end()) {
        throw std::runtime_error("Unknown segment source: " + seg.source);
      }
      sourceMap = &it->second;
    }

    auto sit = sourceMap->find(seg.seqName);
    if (sit == sourceMap->end()) {
      throw std::runtime_error("Sequence not found: " + seg.seqName + " from source " + seg.source);
    }
    std::string seq = sit->second;
    if (seg.reverse) {
      seq = reverseComplement(seq);
    }
    if (seg.start < 0 || seg.end <= seg.start || seg.end > static_cast<int>(seq.size())) {
      throw std::runtime_error("Invalid segment range for " + seg.seqName);
    }
    pieceSeqs.push_back(seq.substr(seg.start, seg.end - seg.start));
  }

  std::string merged;
  for (const auto& p : pieceSeqs) {
    merged += p;
  }

  FastaMap out{{request.outputSeqName.empty() ? "stitched" : request.outputSeqName, merged}};
  writeFasta(request.outputFastaPath, out);

  StitchResult result;
  result.outputFastaPath = request.outputFastaPath;
  result.outputLogPath = request.outputFastaPath + ".session.json";
  result.mergedLength = merged.size();

  for (std::size_t i = 0; i + 1 < pieceSeqs.size(); ++i) {
    BreakpointSummary s;
    s.index = static_cast<int>(i);
    const auto& left = pieceSeqs[i];
    const auto& right = pieceSeqs[i + 1];
    s.leftFlankMatch = compareSuffix(left, right, 50);
    s.rightFlankMatch = comparePrefix(left, right, 50);
    s.preview = previewJunction(left, right, request.contextBp);
    result.breakpoints.push_back(std::move(s));
  }

  std::ofstream log(result.outputLogPath);
  log << "{\n";
  log << "  \"output_fasta\": \"" << toJsonString(result.outputFastaPath) << "\",\n";
  log << "  \"output_name\": \"" << toJsonString(request.outputSeqName.empty() ? "stitched" : request.outputSeqName)
      << "\",\n";
  log << "  \"merged_length\": " << result.mergedLength << ",\n";
  log << "  \"context_bp\": " << request.contextBp << ",\n";
  log << "  \"segments\": [\n";
  for (std::size_t i = 0; i < request.segments.size(); ++i) {
    const auto& seg = request.segments[i];
    log << "    {\"source\":\"" << toJsonString(seg.source)
        << "\",\"name\":\"" << toJsonString(seg.seqName)
        << "\",\"start\":" << seg.start
        << ",\"end\":" << seg.end
        << ",\"reverse\":" << (seg.reverse ? "true" : "false")
        << "}";
    if (i + 1 < request.segments.size()) {
      log << ',';
    }
    log << "\n";
  }
  log << "  ],\n";
  log << "  \"breakpoints\": [\n";
  for (std::size_t i = 0; i < result.breakpoints.size(); ++i) {
    const auto& bp = result.breakpoints[i];
    log << "    {\"index\":" << bp.index
        << ",\"left_flank_match\":" << (bp.leftFlankMatch ? "true" : "false")
        << ",\"right_flank_match\":" << (bp.rightFlankMatch ? "true" : "false")
        << ",\"preview\":\"" << toJsonString(bp.preview) << "\"}";
    if (i + 1 < result.breakpoints.size()) {
      log << ',';
    }
    log << "\n";
  }
  log << "  ]\n";
  log << "}\n";

  return result;
}

}  // namespace gapneedle
