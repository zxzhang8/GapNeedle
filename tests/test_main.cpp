#include "gapneedle/fasta_io.hpp"
#include "gapneedle/guided_stitch_service.hpp"
#include "gapneedle/mapping_service.hpp"
#include "gapneedle/paf.hpp"
#include "gapneedle/facade.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

int main() {
  {
    const std::string seq = "ACGTN";
    const std::string rc = gapneedle::reverseComplement(seq);
    assert(rc == "NACGT");
  }

  {
    std::ofstream paf("/tmp/gapneedle_test.paf");
    paf << "q1\t100\t10\t40\t+\tt1\t120\t20\t50\t25\t30\t60\tcg:Z:30M\n";
    paf.close();

    auto recs = gapneedle::parsePaf("/tmp/gapneedle_test.paf", "t1", "q1");
    assert(recs.size() == 1);
    auto m = gapneedle::mapQueryToTargetDetail(recs[0], 15);
    assert(m.reason == "ok");
    assert(m.tPos.has_value());
    assert(m.tPos.value() == 25);
  }

  {
    std::ofstream paf("/tmp/gapneedle_guided_test.paf");
    paf << "q1\t100\t0\t20\t+\tt1\t200\t0\t20\t20\t20\t60\tcg:Z:20M\n";
    paf << "q1\t100\t20\t50\t+\tt1\t200\t20\t50\t30\t30\t60\tcg:Z:30M\n";
    paf.close();

    gapneedle::GuidedStitchService guided;
    gapneedle::GuidedSeedRequest seedReq;
    seedReq.pafPath = "/tmp/gapneedle_guided_test.paf";
    seedReq.targetSeq = "t1";
    seedReq.querySeq = "q1";
    auto seed = guided.seedCandidates(seedReq);
    assert(!seed.candidates.empty());
    assert(seed.candidates.front().axisStart == 0);

    gapneedle::GuidedStepRequest nextReq;
    nextReq.pafPath = seedReq.pafPath;
    nextReq.targetSeq = seedReq.targetSeq;
    nextReq.querySeq = seedReq.querySeq;
    nextReq.lastAxisEnd = 20;
    nextReq.constraints.minProgressBp = 1;
    nextReq.chosenPath.push_back(seed.candidates.front().segment);
    auto next = guided.nextCandidates(nextReq);
    assert(!next.exhausted);
    assert(!next.candidates.empty());
    assert(next.candidates.front().axisStart >= 20);
  }

  {
    std::ofstream paf("/tmp/gapneedle_guided_clip_suffix_test.paf");
    paf << "q1\t640\t0\t640\t+\tt1\t640\t0\t640\t640\t640\t60\tcg:Z:640M\n";
    paf.close();

    gapneedle::GuidedStitchService guided;
    gapneedle::GuidedStepRequest req;
    req.pafPath = "/tmp/gapneedle_guided_clip_suffix_test.paf";
    req.targetSeq = "t1";
    req.querySeq = "q1";
    req.lastAxisEnd = 314;
    req.constraints.minProgressBp = 1;
    req.chosenPath.push_back(gapneedle::Segment{"q", "q1", 0, 314, false});
    auto next = guided.nextCandidates(req);
    assert(!next.exhausted);
    bool hasT = false;
    bool hasQ = false;
    for (const auto& c : next.candidates) {
      if (c.segment.source == "t" && c.segment.start == 314 && c.segment.end == 640) {
        hasT = true;
      }
      if (c.segment.source == "q" && c.segment.start == 314 && c.segment.end == 640) {
        hasQ = true;
      }
    }
    assert(hasT);
    assert(hasQ);
  }

#if GAPNEEDLE_HAS_MINIMAP2
  {
    const std::string targetPath = "/tmp/gapneedle_mm2_target.fa";
    const std::string queryPath = "/tmp/gapneedle_mm2_query.fa";
    const std::string pafPath = "/tmp/gapneedle_mm2_cache/deep/result/gapneedle_mm2_result.paf";
    const std::string indexDir = "/tmp/gapneedle_mm2_index_cache";
    const std::string seq =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    {
      std::ofstream target(targetPath);
      target << ">t1\n" << seq << '\n';
    }
    {
      std::ofstream query(queryPath);
      query << ">q1\n" << seq << '\n';
    }

    gapneedle::GapNeedleFacade facade;
    gapneedle::AlignmentRequest req;
    req.targetFasta = targetPath;
    req.queryFasta = queryPath;
    req.targetSeq = "t1";
    req.querySeq = "q1";
    req.outputPafPath = pafPath;
    req.preset = "map-ont";
    req.threads = 1;
    req.useIndexCache = true;
    req.indexCacheDir = indexDir;
    req.reuseExisting = false;

    std::filesystem::remove_all("/tmp/gapneedle_mm2_cache");

    auto first = facade.align(req);
    assert(!first.pafPath.empty());
    assert(!first.skipped);
    assert(!first.warnings.empty());
    assert(std::filesystem::exists(pafPath));

    auto recs = gapneedle::parsePaf(pafPath, "t1", "q1");
    assert(!recs.empty());
    bool hasCg = false;
    for (const auto& extra : recs.front().extras) {
      if (extra.rfind("cg:Z:", 0) == 0) {
        hasCg = true;
      }
    }
    assert(hasCg);

    req.outputPafPath = "/tmp/gapneedle_mm2_result_second.paf";
    auto second = facade.align(req);
    assert(!second.warnings.empty());
    bool sawCacheHit = false;
    for (const auto& line : second.warnings) {
      if (line.find("cache hit") != std::string::npos) {
        sawCacheHit = true;
      }
    }
    assert(sawCacheHit);
  }
#endif

  std::cout << "All tests passed\n";
  return 0;
}
