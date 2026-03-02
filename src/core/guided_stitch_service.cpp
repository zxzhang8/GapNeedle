#include "gapneedle/guided_stitch_service.hpp"

#include "gapneedle/paf.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace gapneedle {

namespace {

struct CandidateBuild {
  GuidedCandidate candidate;
};

std::string segKey(const Segment& s) {
  std::ostringstream oss;
  oss << s.source << ":" << s.seqName << ":" << s.start << ":" << s.end << ":" << (s.reverse ? 1 : 0);
  return oss.str();
}

double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

std::string groupByScore(double score) {
  if (score >= 0.70) return "strong";
  if (score >= 0.45) return "acceptable";
  return "risk";
}

std::vector<CandidateBuild> buildFromRecords(const std::vector<AlignmentRecord>& recs) {
  std::vector<CandidateBuild> out;
  out.reserve(recs.size() * 2);
  for (std::size_t i = 0; i < recs.size(); ++i) {
    const auto& r = recs[i];
    if (r.tEnd <= r.tStart || r.qEnd <= r.qStart) {
      continue;
    }

    const std::string recId = "rec#" + std::to_string(i + 1);

    GuidedCandidate tCand;
    tCand.segment.source = "t";
    tCand.segment.seqName = r.tName;
    tCand.segment.start = r.tStart;
    tCand.segment.end = r.tEnd;
    tCand.segment.reverse = false;
    tCand.axisStart = r.tStart;
    tCand.axisEnd = r.tEnd;
    tCand.unclippedAxisEnd = r.tEnd;
    tCand.supportCount = 1;
    tCand.recordId = recId;
    out.push_back(CandidateBuild{tCand});

    GuidedCandidate qCand;
    qCand.segment.source = "q";
    qCand.segment.seqName = r.qName;
    qCand.segment.start = r.qStart;
    qCand.segment.end = r.qEnd;
    qCand.segment.reverse = (r.strand == '-');
    qCand.axisStart = r.tStart;
    qCand.axisEnd = r.tEnd;
    qCand.unclippedAxisEnd = r.tEnd;
    qCand.supportCount = 1;
    qCand.recordId = recId;
    out.push_back(CandidateBuild{qCand});
  }
  return out;
}

std::vector<CandidateBuild> buildSuffixFromRecords(const std::vector<AlignmentRecord>& recs, int lastAxisEnd) {
  std::vector<CandidateBuild> out;
  out.reserve(recs.size() * 2);
  for (std::size_t i = 0; i < recs.size(); ++i) {
    const auto& r = recs[i];
    if (r.tEnd <= r.tStart || r.qEnd <= r.qStart) {
      continue;
    }
    if (lastAxisEnd >= r.tEnd) {
      continue;
    }

    const int suffixAxisStart = std::max(r.tStart, lastAxisEnd);
    const int suffixAxisEnd = r.tEnd;
    if (suffixAxisEnd <= suffixAxisStart) {
      continue;
    }

    const std::string recId = "rec#" + std::to_string(i + 1);

    GuidedCandidate tCand;
    tCand.segment.source = "t";
    tCand.segment.seqName = r.tName;
    tCand.segment.start = suffixAxisStart;
    tCand.segment.end = suffixAxisEnd;
    tCand.segment.reverse = false;
    tCand.axisStart = suffixAxisStart;
    tCand.axisEnd = suffixAxisEnd;
    tCand.unclippedAxisEnd = r.tEnd;
    tCand.supportCount = 1;
    tCand.recordId = recId;
    tCand.rationale = "suffix candidate from target axis";
    out.push_back(CandidateBuild{tCand});

    const int tSpan = r.tEnd - r.tStart;
    const int qSpan = r.qEnd - r.qStart;
    if (tSpan <= 0 || qSpan <= 0) {
      continue;
    }
    auto mapAxisToQuery = [&](int axisPos) {
      const double ratio = static_cast<double>(axisPos - r.tStart) / static_cast<double>(tSpan);
      int q = r.qStart + static_cast<int>(std::lround(ratio * static_cast<double>(qSpan)));
      q = std::max(r.qStart, std::min(r.qEnd, q));
      return q;
    };
    int qStart = mapAxisToQuery(suffixAxisStart);
    int qEnd = mapAxisToQuery(suffixAxisEnd);
    if (qEnd <= qStart) {
      qEnd = std::min(r.qEnd, qStart + 1);
    }
    if (qEnd <= qStart) {
      continue;
    }

    GuidedCandidate qCand;
    qCand.segment.source = "q";
    qCand.segment.seqName = r.qName;
    qCand.segment.start = qStart;
    qCand.segment.end = qEnd;
    qCand.segment.reverse = (r.strand == '-');
    qCand.axisStart = suffixAxisStart;
    qCand.axisEnd = suffixAxisEnd;
    qCand.unclippedAxisEnd = r.tEnd;
    qCand.supportCount = 1;
    qCand.recordId = recId;
    qCand.rationale = "suffix candidate mapped from target axis";
    out.push_back(CandidateBuild{qCand});
  }
  return out;
}

void mergeDuplicates(std::vector<CandidateBuild>* inout) {
  std::unordered_map<std::string, CandidateBuild> merged;
  for (const auto& b : *inout) {
    const std::string key = segKey(b.candidate.segment) + "|" + std::to_string(b.candidate.axisStart) + ":" +
                            std::to_string(b.candidate.axisEnd);
    auto it = merged.find(key);
    if (it == merged.end()) {
      merged.emplace(key, b);
      continue;
    }
    it->second.candidate.supportCount += 1;
    if (!it->second.candidate.recordId.empty()) {
      it->second.candidate.recordId += "," + b.candidate.recordId;
    } else {
      it->second.candidate.recordId = b.candidate.recordId;
    }
  }
  inout->clear();
  inout->reserve(merged.size());
  for (auto& [_, b] : merged) {
    inout->push_back(std::move(b));
  }
}

double seedScore(const GuidedCandidate& c, const GuidedConstraints& cfg) {
  const double d = static_cast<double>(std::max(0, c.axisStart));
  const double near = 1.0 - clamp01(d / std::max(1.0, static_cast<double>(cfg.nearZeroWindow)));
  const double support = clamp01(static_cast<double>(c.supportCount) / 4.0);
  const double len = clamp01(static_cast<double>(std::max(0, c.axisEnd - c.axisStart)) / 200000.0);
  return 0.45 * near + 0.35 * support + 0.20 * len;
}

double stepScore(const GuidedCandidate& c, int lastAxisEnd, const GuidedConstraints& cfg) {
  const int progressBp = std::max(0, c.axisEnd - lastAxisEnd);
  const int jumpBp = std::max(0, c.axisStart - lastAxisEnd);
  const double progress = clamp01(static_cast<double>(progressBp) / 300000.0);
  const double jumpPenalty = clamp01(static_cast<double>(jumpBp) / std::max(1.0, static_cast<double>(cfg.maxJumpBp)));
  const double support = clamp01(static_cast<double>(c.supportCount) / 4.0);
  const double len = clamp01(static_cast<double>(std::max(0, c.axisEnd - c.axisStart)) / 200000.0);
  return 0.40 * progress + 0.25 * support + 0.20 * len + 0.15 * (1.0 - jumpPenalty);
}

template <typename T>
void sortAndTrim(std::vector<T>* items, int topK) {
  std::sort(items->begin(), items->end(), [](const T& a, const T& b) {
    if (std::abs(a.candidate.score - b.candidate.score) > 1e-12) {
      return a.candidate.score > b.candidate.score;
    }
    if (a.candidate.axisStart != b.candidate.axisStart) {
      return a.candidate.axisStart < b.candidate.axisStart;
    }
    return a.candidate.axisEnd > b.candidate.axisEnd;
  });
  if (topK > 0 && static_cast<int>(items->size()) > topK) {
    items->resize(static_cast<std::size_t>(topK));
  }
}

}  // namespace

GuidedSeedResult GuidedStitchService::seedCandidates(const GuidedSeedRequest& request) const {
  if (request.pafPath.empty() || request.targetSeq.empty() || request.querySeq.empty()) {
    throw std::runtime_error("guided seed requires --paf --target-seq --query-seq");
  }

  auto recs = parsePaf(request.pafPath, request.targetSeq, request.querySeq);
  if (recs.empty()) {
    return GuidedSeedResult{{}, {"no records found in PAF for selected sequences"}};
  }

  auto built = buildFromRecords(recs);
  mergeDuplicates(&built);

  std::vector<CandidateBuild> strict;
  std::vector<CandidateBuild> fallback;
  for (auto& b : built) {
    if (b.candidate.axisStart == 0) {
      strict.push_back(b);
    } else if (b.candidate.axisStart > 0 && b.candidate.axisStart <= request.constraints.nearZeroWindow) {
      b.candidate.fallbackNearZero = true;
      fallback.push_back(b);
    }
  }

  std::vector<std::string> warnings;
  auto* active = &strict;
  if (strict.empty()) {
    active = &fallback;
    if (!fallback.empty()) {
      warnings.emplace_back("no strict seed at axis 0, fallback to near-zero candidates");
    } else {
      warnings.emplace_back("no seed candidates found at or near axis 0");
    }
  }

  for (auto& b : *active) {
    b.candidate.score = seedScore(b.candidate, request.constraints);
    b.candidate.group = groupByScore(b.candidate.score);
    b.candidate.rationale = b.candidate.fallbackNearZero
                                ? "near-zero start fallback with alignment support"
                                : "strict axis-0 start candidate";
  }

  sortAndTrim(active, request.maxSeeds);

  GuidedSeedResult out;
  out.warnings = std::move(warnings);
  out.candidates.reserve(active->size());
  for (auto& b : *active) {
    out.candidates.push_back(std::move(b.candidate));
  }
  return out;
}

GuidedStepResult GuidedStitchService::nextCandidates(const GuidedStepRequest& request) const {
  if (request.pafPath.empty() || request.targetSeq.empty() || request.querySeq.empty()) {
    throw std::runtime_error("guided next requires --paf --target-seq --query-seq");
  }

  auto recs = parsePaf(request.pafPath, request.targetSeq, request.querySeq);
  if (recs.empty()) {
    return GuidedStepResult{{}, true, {"no records found in PAF for selected sequences"}};
  }

  auto built = buildSuffixFromRecords(recs, request.lastAxisEnd);
  mergeDuplicates(&built);

  std::unordered_set<std::string> used;
  for (const auto& s : request.chosenPath) {
    used.insert(segKey(s));
  }
  Segment lastSelected;
  bool hasLastSelected = false;
  if (!request.chosenPath.empty()) {
    lastSelected = request.chosenPath.back();
    hasLastSelected = true;
  }

  std::vector<CandidateBuild> next;
  next.reserve(built.size());
  for (auto& b : built) {
    const auto& c = b.candidate;
    if (used.count(segKey(c.segment)) > 0) {
      continue;
    }
    if (hasLastSelected && c.segment.source == lastSelected.source && c.segment.seqName == lastSelected.seqName &&
        c.segment.end <= lastSelected.end) {
      continue;
    }
    if (c.axisStart < request.lastAxisEnd) {
      continue;  // strict monotonic progression
    }
    const int jumpBp = c.axisStart - request.lastAxisEnd;
    if (jumpBp > request.constraints.maxJumpBp) {
      continue;
    }
    const int progressBp = c.axisEnd - request.lastAxisEnd;
    if (progressBp < request.constraints.minProgressBp) {
      continue;
    }
    b.candidate.score = stepScore(c, request.lastAxisEnd, request.constraints);
    b.candidate.group = groupByScore(b.candidate.score);
    std::ostringstream reason;
    reason << "progress " << progressBp << "bp, jump " << jumpBp << "bp, support " << c.supportCount;
    b.candidate.rationale = reason.str();
    next.push_back(std::move(b));
  }

  sortAndTrim(&next, request.maxNext);

  GuidedStepResult out;
  out.exhausted = next.empty();
  if (out.exhausted) {
    out.warnings.emplace_back("no monotonic next candidate found");
  }
  out.candidates.reserve(next.size());
  for (auto& b : next) {
    out.candidates.push_back(std::move(b.candidate));
  }
  return out;
}

}  // namespace gapneedle
