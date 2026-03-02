#include "gapneedle/facade.hpp"
#include "gapneedle/telomere_service.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using gapneedle::AlignmentRequest;
using gapneedle::GapNeedleFacade;
using gapneedle::Segment;
using gapneedle::StitchRequest;

namespace {

std::unordered_map<std::string, std::vector<std::string>> parseArgs(int argc, char* argv[]) {
  std::unordered_map<std::string, std::vector<std::string>> opts;
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    if (key.rfind("--", 0) != 0) {
      continue;
    }
    if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
      opts[key].push_back(argv[++i]);
    } else {
      opts[key].push_back("true");
    }
  }
  return opts;
}

std::string getOne(const std::unordered_map<std::string, std::vector<std::string>>& opts,
                   const std::string& key,
                   const std::string& def = "") {
  auto it = opts.find(key);
  if (it == opts.end() || it->second.empty()) {
    return def;
  }
  return it->second.back();
}

std::vector<std::string> getMany(const std::unordered_map<std::string, std::vector<std::string>>& opts,
                                 const std::string& key) {
  auto it = opts.find(key);
  if (it == opts.end()) {
    return {};
  }
  return it->second;
}

void printUsage() {
  std::cout << "gapneedle_cli --cmd <align|stitch|scan-gaps|check-telomere|guided-seed|guided-next> [options]\n"
            << "  align: --target-fasta --query-fasta --target-seq --query-seq [--output] [--preset] [--threads]\n"
            << "  stitch: --target-fasta --query-fasta --output --segment src:name:start:end[:rc] (repeatable)\n"
            << "  scan-gaps: --target-fasta [--min-gap]\n"
            << "  check-telomere: --target-fasta --seq-name\n"
            << "  guided-seed: --paf --target-seq --query-seq [--max-seeds] [--near-zero-window]\n"
            << "  guided-next: --paf --target-seq --query-seq --last-axis-end [--max-next] [--max-jump-bp] [--min-progress-bp]\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  auto opts = parseArgs(argc, argv);
  const std::string cmd = getOne(opts, "--cmd");
  if (cmd.empty() || getOne(opts, "--help") == "true") {
    printUsage();
    return cmd.empty() ? 2 : 0;
  }

  GapNeedleFacade facade;
  try {
    if (cmd == "align") {
      AlignmentRequest req;
      req.targetFasta = getOne(opts, "--target-fasta");
      req.queryFasta = getOne(opts, "--query-fasta");
      req.targetSeq = getOne(opts, "--target-seq");
      req.querySeq = getOne(opts, "--query-seq");
      req.outputPafPath = getOne(opts, "--output");
      req.preset = getOne(opts, "--preset", "asm10");
      req.threads = std::stoi(getOne(opts, "--threads", "4"));
      auto r = facade.align(req);
      std::cout << "PAF: " << r.pafPath << "\n";
      std::cout << "Skipped: " << (r.skipped ? "true" : "false") << "\n";
    } else if (cmd == "stitch") {
      StitchRequest req;
      req.targetFasta = getOne(opts, "--target-fasta");
      req.queryFasta = getOne(opts, "--query-fasta");
      req.outputFastaPath = getOne(opts, "--output");
      req.outputSeqName = getOne(opts, "--output-name", "stitched");

      for (const auto& s : getMany(opts, "--segment")) {
        std::vector<std::string> parts;
        std::size_t start = 0;
        while (start <= s.size()) {
          auto p = s.find(':', start);
          if (p == std::string::npos) {
            parts.push_back(s.substr(start));
            break;
          }
          parts.push_back(s.substr(start, p - start));
          start = p + 1;
        }
        if (parts.size() < 4) {
          continue;
        }
        Segment seg;
        seg.source = parts[0];
        seg.seqName = parts[1];
        seg.start = std::stoi(parts[2]);
        seg.end = std::stoi(parts[3]);
        seg.reverse = (parts.size() > 4 && parts[4] == "rc");
        req.segments.push_back(seg);
      }

      auto r = facade.stitch(req);
      std::cout << "Output FASTA: " << r.outputFastaPath << "\n";
      std::cout << "Session log: " << r.outputLogPath << "\n";
      std::cout << "Merged length: " << r.mergedLength << "\n";
    } else if (cmd == "scan-gaps") {
      auto gaps = facade.scanGaps(getOne(opts, "--target-fasta"), std::stoi(getOne(opts, "--min-gap", "10")));
      for (const auto& [name, s, e] : gaps) {
        std::cout << name << "\t" << s << "\t" << e << "\n";
      }
    } else if (cmd == "check-telomere") {
      auto [left, right] = gapneedle::checkTelomere(getOne(opts, "--target-fasta"), getOne(opts, "--seq-name"));
      std::cout << "left=" << (left ? "true" : "false") << " right=" << (right ? "true" : "false") << "\n";
    } else if (cmd == "guided-seed") {
      gapneedle::GuidedSeedRequest req;
      req.pafPath = getOne(opts, "--paf");
      req.targetSeq = getOne(opts, "--target-seq");
      req.querySeq = getOne(opts, "--query-seq");
      req.maxSeeds = std::stoi(getOne(opts, "--max-seeds", "12"));
      req.constraints.nearZeroWindow = std::stoi(getOne(opts, "--near-zero-window", "1000"));

      auto r = facade.guidedSeed(req);
      for (const auto& w : r.warnings) {
        std::cout << "warning\t" << w << "\n";
      }
      for (std::size_t i = 0; i < r.candidates.size(); ++i) {
        const auto& c = r.candidates[i];
        std::cout << i
                  << "\t" << c.segment.source << ":" << c.segment.seqName << ":" << c.segment.start << ":" << c.segment.end
                  << (c.segment.reverse ? ":rc" : "")
                  << "\taxis=" << c.axisStart << "-" << c.axisEnd
                  << "\tscore=" << c.score
                  << "\tgroup=" << c.group
                  << "\tsupport=" << c.supportCount
                  << "\t" << c.rationale << "\n";
      }
    } else if (cmd == "guided-next") {
      gapneedle::GuidedStepRequest req;
      req.pafPath = getOne(opts, "--paf");
      req.targetSeq = getOne(opts, "--target-seq");
      req.querySeq = getOne(opts, "--query-seq");
      req.lastAxisEnd = std::stoi(getOne(opts, "--last-axis-end", "0"));
      req.maxNext = std::stoi(getOne(opts, "--max-next", "12"));
      req.constraints.maxJumpBp = std::stoi(getOne(opts, "--max-jump-bp", "200000"));
      req.constraints.minProgressBp = std::stoi(getOne(opts, "--min-progress-bp", "200"));

      auto r = facade.guidedNext(req);
      std::cout << "exhausted=" << (r.exhausted ? "true" : "false") << "\n";
      for (const auto& w : r.warnings) {
        std::cout << "warning\t" << w << "\n";
      }
      for (std::size_t i = 0; i < r.candidates.size(); ++i) {
        const auto& c = r.candidates[i];
        std::cout << i
                  << "\t" << c.segment.source << ":" << c.segment.seqName << ":" << c.segment.start << ":" << c.segment.end
                  << (c.segment.reverse ? ":rc" : "")
                  << "\taxis=" << c.axisStart << "-" << c.axisEnd
                  << "\tscore=" << c.score
                  << "\tgroup=" << c.group
                  << "\tsupport=" << c.supportCount
                  << "\t" << c.rationale << "\n";
      }
    } else {
      printUsage();
      return 2;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
