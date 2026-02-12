#include "gapneedle/mapping_service.hpp"

#include <cctype>
#include <vector>

namespace gapneedle {

namespace {

std::string getTagValue(const AlignmentRecord& rec, const std::string& tag) {
  const std::string prefix = tag + ":";
  for (const auto& e : rec.extras) {
    if (e.rfind(prefix, 0) == 0) {
      auto p1 = e.find(':');
      if (p1 == std::string::npos) {
        continue;
      }
      auto p2 = e.find(':', p1 + 1);
      if (p2 == std::string::npos) {
        continue;
      }
      return e.substr(p2 + 1);
    }
  }
  return {};
}

std::vector<std::pair<int, char>> parseCigar(const std::string& cigar) {
  std::vector<std::pair<int, char>> ops;
  int num = 0;
  bool hasNum = false;
  for (char ch : cigar) {
    if (std::isdigit(static_cast<unsigned char>(ch))) {
      hasNum = true;
      num = num * 10 + (ch - '0');
      continue;
    }
    if (hasNum) {
      ops.emplace_back(num, ch);
      num = 0;
      hasNum = false;
    }
  }
  return ops;
}

}  // namespace

MappingResult mapQueryToTargetDetail(const AlignmentRecord& rec, int qPos) {
  MappingResult result;
  result.reason = "no_mapping";
  result.qPos = qPos;

  for (char op : std::string("M=XIDNSHP")) {
    result.countsBefore[op] = 0;
    result.countsTotal[op] = 0;
  }

  const std::string cigar = getTagValue(rec, "cg");
  if (cigar.empty()) {
    result.reason = "missing_cigar";
    return result;
  }
  if (qPos < rec.qStart || qPos >= rec.qEnd) {
    result.reason = "out_of_range";
    return result;
  }

  int qPosOriented = qPos;
  int qCursor = rec.qStart;
  if (rec.strand == '-') {
    qPosOriented = rec.qLen - 1 - qPos;
    qCursor = rec.qLen - rec.qEnd;
  }
  result.qPosOriented = qPosOriented;
  int tCursor = rec.tStart;

  const auto ops = parseCigar(cigar);
  for (const auto& [len, op] : ops) {
    if (result.countsTotal.count(op)) {
      result.countsTotal[op] += len;
    }
    if (op == 'M' || op == '=' || op == 'X') {
      if (qPosOriented < qCursor + len) {
        result.tPos = tCursor + (qPosOriented - qCursor);
        result.reason = "ok";
        result.op = op;
        result.opLen = len;
        result.opOffset = qPosOriented - qCursor;
        return result;
      }
      qCursor += len;
      tCursor += len;
      result.countsBefore[op] += len;
      result.qConsumedBefore += len;
      result.tConsumedBefore += len;
    } else if (op == 'I' || op == 'S') {
      if (qPosOriented < qCursor + len) {
        result.reason = "insertion";
        result.op = op;
        result.opLen = len;
        result.opOffset = qPosOriented - qCursor;
        return result;
      }
      qCursor += len;
      result.countsBefore[op] += len;
      result.qConsumedBefore += len;
    } else if (op == 'D' || op == 'N') {
      tCursor += len;
      result.countsBefore[op] += len;
      result.tConsumedBefore += len;
    } else if (op == 'H' || op == 'P') {
      continue;
    } else {
      result.reason = "bad_cigar";
      result.op = op;
      result.opLen = len;
      return result;
    }
  }

  result.reason = "no_mapping";
  return result;
}

}  // namespace gapneedle
